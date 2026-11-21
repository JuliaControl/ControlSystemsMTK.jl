using ControlSystemsMTK, ControlSystemsBase, ModelingToolkit, RobustAndOptimalControl, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D, connect
using ControlSystemsMTK: lifted_matrices

# A continuous-time first-order plant, a filter on a fast clock and a PI controller on a slow
# clock. The periods are exactly representable such that the ticks of the clocks coincide.
dtf = 0.125
ratio = 3
dts = ratio * dtf
kf = ShiftIndex(Clock(dtf))
ks = ShiftIndex(Clock(dts))
kpv, kiv, gv = 1.5, 1.0, 1.0
@variables x(t) y(t) d(t) r(t) e(t) xi(t) ud(t) w(t)
@variables uh(t) [output = true] u(t) [input = true]
@variables wo(t) [output = true] wi(t) [input = true] w2(t)
@parameters kp = kpv ki = kiv g = gv

plant = [
    D(x) ~ -x + u + d
    y ~ x
]
fastfilter = [w(kf) ~ 0.5 * w(kf - 1) + 0.5 * Sample(dtf)(y)]
# `measurement` is a signal on the slow clock
slow_controller(measurement) = [
    e ~ r - measurement
    xi(ks) ~ xi(ks - 1) + ki * dts * e(ks)
    ud ~ kp * e + xi
]
hold = [uh ~ Hold(ud); connect(uh, :plant_input, u)]
op = Dict(x => 0.0, d => 0.0, r => 0.0, xi => 0.0, w => 0.0)

slow(signal) = Sample(dts)(Hold(signal))
@named loop = System([plant; fastfilter; slow_controller(slow(w)); hold], t)
# The gain `g` and an analysis point in the fast partition
@named loop_fast_break = System([plant; fastfilter; wo ~ g * w; connect(wo, :fast_signal, wi); w2(kf) ~ wi(kf); slow_controller(slow(w2)); hold; r ~ 0; d ~ 0], t)

# Reference simulation of `loop`, written directly from the equations. `dfun(m)` is the
# value of `d` that is held from fast tick `m`, `rfun(m)` the value of `r` at slow tick `m`.
# Returns the outputs `y`, `w` at the fast ticks and `ud` at the slow ticks.
function simulate_reference(dfun, rfun, nperiods; gain = gv, x0 = zeros(5))
    a = exp(-dtf)
    xv, wprev, xiprev, uheld = x0[1], x0[2], x0[3], x0[4]
    ys, ws, uds = zeros(ComplexF64, 0), zeros(ComplexF64, 0), zeros(ComplexF64, 0)
    for m in 0:ratio*nperiods-1
        yv = xv
        wv = 0.5 * wprev + 0.5 * yv
        if m % ratio == 0
            ev = rfun(m ÷ ratio) - gain * wv
            xiv = xiprev + kiv * dts * ev
            uheld = kpv * ev + xiv
            xiprev = xiv
            push!(uds, uheld)
        end
        wprev = wv
        push!(ys, yv); push!(ws, wv)
        xv = a * xv + (1 - a) * (uheld + dfun(m))
    end
    ys, ws, uds
end

if isdefined(ModelingToolkit, :linearize_hybrid)

ml = multirate_linearization(loop, [d, r], [y, w, ud]; op)

@testset "Structure" begin
    @test ml isa MultirateLinearization
    @test ml.Tc ≈ dts
    @test ml.event_times ≈ dtf .* (0:ratio-1)
    @test [length(p.times) for p in ml.inputs] == [ratio, 1]
    @test [length(p.times) for p in ml.outputs] == [ratio, ratio, 1]
    G = lifted_ss(ml)
    @test G.Ts ≈ dts
    @test G.nu == ratio + 1
    @test G.ny == 2ratio + 1
    @test G.u[1] == Symbol("d(t)@1")
    @test G.u[end] == Symbol("r(t)")
    # The feedthrough matrix of a causal model is lower triangular in chronological order
    F, Gm, H, E = lifted_matrices(ml)
    yd = [p.channels for p in ml.outputs]
    @test iszero(E[yd[1][1], ml.inputs[1].channels[2:end]])
    @test iszero(E[yd[2][2], ml.inputs[1].channels[3:end]])
end

@testset "Lifted system and reference simulation" begin
    nper = 6
    dv = randn(ratio * nper)
    rv = randn(nper)
    ys, ws, uds = simulate_reference(m -> dv[m+1], m -> rv[m+1], nper)
    G = lifted_ss(ml)
    U = [reshape(dv, ratio, nper); rv']
    res = lsim(G.sys, U)
    @test vec(res.y[ml.outputs[1].channels, :]) ≈ real(ys)
    @test vec(res.y[ml.outputs[2].channels, :]) ≈ real(ws)
    @test vec(res.y[ml.outputs[3].channels, :]) ≈ real(uds)
end

@testset "Multipliers" begin
    λ = multipliers(ml)
    @test maximum(abs, λ) < 1
    # The multipliers do not depend on the start of the period
    ml2 = multirate_linearization(loop, [d, r], [y, w, ud]; op, t0 = dtf)
    @test sort(abs.(multipliers(ml2))) ≈ sort(abs.(λ)) atol = 1e-8
    @test ml2.event_times ≈ dtf .* (1:ratio)
end

@testset "Harmonic response" begin
    wv = [0.3, 2.0, 7.0, 20.0]
    Ω = 2π / dts
    nper = 400
    for (i, w0) in enumerate(wv)
        # Steady-state response to d = exp(im*w0*t), held on the fast grid
        ys, ws, uds = simulate_reference(m -> cis(w0 * m * dtf), m -> 0, nper)
        tail = (nper - 50) * ratio + 1 : nper * ratio
        tf_ = dtf .* (tail .- 1)
        for k in 0:ratio-1
            c = sum(ys[tail] .* cis.(-(w0 + k * Ω) .* tf_)) / length(tail)
            @test harmonic_response(ml, "d(t)", "y(t)", [w0], k)[] ≈ c atol = 1e-6
        end
        @test fundamental_response(ml, 1, 1, [w0])[] ≈ harmonic_response(ml, 1, 1, [w0], 0)[]
        # The output on the slow clock
        tails = nper - 50 + 1 : nper
        c = sum(uds[tails] .* cis.(-w0 .* dts .* (tails .- 1))) / length(tails)
        @test fundamental_response(ml, "d(t)", "ud(t)", [w0])[] ≈ c atol = 1e-6
        # Parseval: the performance frequency gain collects all harmonics
        pfg = performance_frequency_gain(ml, 1, 1, [w0])[]
        @test pfg^2 ≈ sum(abs2(harmonic_response(ml, 1, 1, [w0], k)[]) for k in 0:ratio-1)
        @test abs(fundamental_response(ml, 1, 1, [w0])[]) <= pfg + 1e-12
        @test pfg <= robust_frequency_gain(ml, 1, 1, [w0])[] + 1e-12
    end
    @test robust_frequency_gain(ml, 1, 1, [0.3])[] ≈ robust_frequency_gain(ml, 1, 1, [0.3 + Ω])[]
    @test_throws "No input named" fundamental_response(ml, "nonexistent", 1, wv)
end

@testset "A model with one clock" begin
    @named single = System([plant; slow_controller(Sample(dts)(y)); hold; r ~ 0], t)
    op1 = Dict(x => 0.0, d => 0.0, xi => 0.0)
    G = named_ss(single, [d], [y]; op = op1, hybrid = true)
    ml1 = multirate_linearization(single, [d], [y]; op = op1)
    @test length(ml1.event_times) == 1
    wv = exp10.(range(-2, log10(π / dts), length = 30))
    @test freqresp(lifted_ss(ml1).sys, wv) ≈ freqresp(G.sys, wv)
    @test fundamental_response(ml1, 1, 1, wv) ≈ vec(freqresp(G.sys, wv))
    @test sort(filter(λ -> abs(λ) > 1e-10, multipliers(ml1)), by = abs) ≈ sort(filter(λ -> abs(λ) > 1e-10, poles(G.sys)), by = abs)
    # Oversampling exposes the behavior between the ticks and leaves the dynamics unchanged
    ml4 = multirate_linearization(single, [d], [y]; op = op1, oversample = 4)
    @test length(ml4.inputs[1].times) == 4
    @test maximum(abs, multipliers(ml4)) ≈ maximum(abs, multipliers(ml1))
    # `single_rate_ss` at the rate of the clock is the same system
    Gs = single_rate_ss(single, [d], [y]; Ts = dts, method = :zoh, op = op1)
    @test freqresp(Gs.sys, wv) ≈ freqresp(G.sys, wv)
    # In continuous time, the controller is converted with the Tustin approximation
    Gc = single_rate_ss(single, [d], [y]; Ts = 0, op = op1)
    z = tf("z", dts)
    Cc = d2c(ss(kpv + kiv * dts * z / (z - 1)), :tustin)
    P = ss(-1.0, 1.0, 1.0, 0.0)
    @test freqresp(Gc.sys, wv) ≈ freqresp(feedback(P, Cc), wv)
    # A conversion that introduces feedthrough in the plant is not rejected
    Gt = single_rate_ss(single, [d], [y]; Ts = dts / 2, op = op1)
    @test Gt.Ts ≈ dts / 2
    @test maximum(abs, poles(Gt.sys)) < 1
end

@testset "Clocked partitions that read each other at the same instant" begin
    # Neither output depends directly on the other, so the model is causal
    @variables p(t) q(t)
    eqs = [
        D(x) ~ -x + d
        p(kf) ~ 0.5 * p(kf - 1) + 0.2 * Sample(dtf)(Hold(q))
        q(ks) ~ 0.4 * q(ks - 1) + 0.3 * Sample(dts)(Hold(p))
    ]
    @named mutual = System(eqs, t)
    @test_throws "not causal" multirate_linearization(mutual, [d], [x]; op = Dict(x => 0.0, d => 0.0, p => 0.0, q => 0.0))
    # The fast partition uses the value of `q` at its previous tick
    @variables qs(t)
    eqs = [
        D(x) ~ -x + d
        qs ~ Sample(dtf)(Hold(q))
        p(kf) ~ 0.5 * p(kf - 1) + 0.2 * qs(kf - 1)
        q(ks) ~ 0.4 * q(ks - 1) + 0.3 * Sample(dts)(Hold(p))
    ]
    @named mutual = System(eqs, t)
    mlm = multirate_linearization(mutual, [d], [x]; op = Dict(x => 0.0, d => 0.0, p => 0.0, q => 0.0, qs => 0.0))
    # One period of the clocked equations, with the state (p, qs, q) of the previous ticks
    function period(v)
        pprev, qsprev, qprev = v
        for m in 0:ratio-1
            pv = 0.5 * pprev + 0.2 * qsprev
            m == 0 && (qprev = 0.4 * qprev + 0.3 * pv)
            pprev, qsprev = pv, qprev
        end
        [pprev, qsprev, qprev]
    end
    Mref = reduce(hcat, period.(eachcol(Matrix(1.0I, 3, 3))))
    nonzero(λ) = sort(filter(l -> abs(l) > 1e-8, λ), by = abs)
    # The plant is not connected to the clocked partitions and contributes its own multiplier
    @test nonzero(multipliers(mlm)) ≈ nonzero([eigvals(Mref); exp(-dts)])
end

@testset "Held signal with direct feedthrough to a sampled signal" begin
    # ud[k] = -a*(x[k] + ud[k-1]): the sampled signal contains the previously held value
    a = 0.5
    eqs = [
        D(x) ~ -x + Hold(ud) + d
        y ~ x + Hold(ud)
        ud ~ -a * Sample(dts)(y)
    ]
    @named feedthrough = System(eqs, t)
    mlf = multirate_linearization(feedthrough, [d], [x]; op = Dict(x => 0.0, d => 0.0))
    b = exp(-dts)
    # State: x and the held ud
    Aref = [b-(1-b)*a -(1-b)*a; -a -a]
    @test sort(filter(λ -> abs(λ) > 1e-10, multipliers(mlf)), by = real) ≈ sort(eigvals(Aref), by = real)
end

@testset "Margins" begin
    @testset "Break point with one instant per period" begin
        @named single_break = System([plant; wo ~ Sample(dts)(y); connect(wo, :measurement, wi); slow_controller(wi); hold; r ~ 0; d ~ 0], t)
        opb = Dict(x => 0.0, xi => 0.0)
        m = multirate_margins(single_break, single_break.measurement; op = opb, N = 2000)
        z = tf("z", dts)
        L = c2d(ss(-1.0, 1.0, 1.0, 0.0), dts) * ss(kpv + kiv * dts * z / (z - 1))
        dm = diskmargin(L)
        @test m.disk_margin_lti_lower_bound.α ≈ dm.α rtol = 1e-4
        @test m.disk_margin_dynamic.α ≈ dm.α rtol = 1e-4
        # The upper limit of the gain is where a closed-loop pole reaches the unit circle
        gmax = m.gain_interval[2]
        @test maximum(abs, poles(feedback(gmax * L))) ≈ 1 atol = 1e-6
        @test m.gain_interval[1] == 0
    end
    @testset "Break point in the fast partition" begin
        opb = Dict(x => 0.0, xi => 0.0, w => 0.0)
        m = multirate_margins(loop_fast_break, loop_fast_break.fast_signal; op = opb)
        gmin, gmax = m.gain_interval
        @test gmin == 0
        @test 1 < gmax < Inf
        # The model with the gain at the limit has a multiplier on the unit circle
        mlg = multirate_linearization(loop_fast_break, [], [x]; op = merge(opb, Dict(g => gmax)))
        @test maximum(abs, multipliers(mlg)) ≈ 1 atol = 1e-6
        mlg = multirate_linearization(loop_fast_break, [], [x]; op = merge(opb, Dict(g => 0.98gmax)))
        @test maximum(abs, multipliers(mlg)) < 1
        # Larger classes of perturbations give smaller margins
        @test m.disk_margin_dynamic.α <= m.disk_margin_lti_lower_bound.α + 1e-6
        @test m.disk_margin_lti_lower_bound.γmax <= gmax + 1e-6
    end
end

@testset "Unsupported models" begin
    @test_throws "multirate_linearization" named_ss(loop, [d, r], [y]; op, hybrid = true)
    kn = ShiftIndex(Clock(dtf * sqrt(2)))
    @variables v(t)
    eqs = [
        D(x) ~ -x + Hold(w) + Hold(v) + d
        w(kf) ~ 0.5 * w(kf - 1) - Sample(dtf)(x)
        v(kn) ~ 0.5 * v(kn - 1) - Sample(dtf * sqrt(2))(x)
    ]
    @named incommensurate = System(eqs, t)
    opi = Dict(x => 0.0, d => 0.0, w => 0.0, v => 0.0)
    @test_throws "not commensurate" multirate_linearization(incommensurate, [d], [x]; op = opi)
    # The approximation by rate conversion is available
    Gc = single_rate_ss(incommensurate, [d], [x]; Ts = 0, op = opi)
    @test Gc.nx == 3
    @test ControlSystemsBase.iscontinuous(Gc)
end

end
