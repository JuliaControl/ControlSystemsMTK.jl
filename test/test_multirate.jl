using ControlSystemsMTK, ControlSystemsBase, ModelingToolkit, RobustAndOptimalControl, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D, connect
using ControlSystemsMTK: lifted_matrices

# A continuous-time first-order plant, a filter on a fast clock and a PI controller on a slow
# clock, which reads the filter through `Sample(Hold(w))`. The periods are exactly
# representable such that the ticks of the clocks coincide.
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
            # `Sample(Hold(w))` has the value of w that was held before this tick
            ev = rfun(m ÷ ratio) - gain * wprev
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

# The tests below require `ModelingToolkit.linearize_hybrid`
@test isdefined(ModelingToolkit, :linearize_hybrid)
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
    @test isstable(ml)
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
    # Each partition reads the other through `Sample(Hold(x))` and thus obtains the value
    # that was held before the tick, also when both tick at the same instant
    @variables p(t) q(t)
    eqs = [
        D(x) ~ -x + d
        p(kf) ~ 0.5 * p(kf - 1) + 0.2 * Sample(dtf)(Hold(q))
        q(ks) ~ 0.4 * q(ks - 1) + 0.3 * Sample(dts)(Hold(p))
    ]
    @named mutual = System(eqs, t)
    mlm = multirate_linearization(mutual, [d], [x]; op = Dict(x => 0.0, d => 0.0, p => 0.0, q => 0.0))
    # One period of the clocked equations, with the state (p, q) of the previous ticks
    function period(v)
        pv, qv = v
        for m in 0:ratio-1
            pnew = 0.5 * pv + 0.2 * qv
            m == 0 && (qv = 0.4 * qv + 0.3 * pv)
            pv = pnew
        end
        [pv, qv]
    end
    Mref = reduce(hcat, period.(eachcol(Matrix(1.0I, 2, 2))))
    nonzero(λ) = sort(filter(l -> abs(l) > 1e-8, λ), by = abs)
    # The plant is not connected to the clocked partitions and contributes its own multiplier
    @test nonzero(multipliers(mlm)) ≈ nonzero([eigvals(Mref); exp(-dts)])

    # Outputs that depend directly on each other at the same instant are not causal. Such a
    # dependency arises with operators that read the value computed at the instant.
    systems = Dict(1 => ss(0.5, 1, 1, 1, dtf), 2 => ss(0.5, 1, 1, 1, dtf))
    source = Dict((1, 1) => (2, 1), (2, 1) => (1, 1))
    @test_throws "not causal" ControlSystemsMTK.clocked_outputs([1, 2], systems, source, (i, j) -> nothing, i -> zeros(1, 2), 2)
    # A dependency in one direction is resolved
    systems = Dict(1 => ss(0.5, 1, 1, 0, dtf), 2 => ss(0.5, 1, 1, 1, dtf))
    ynew = ControlSystemsMTK.clocked_outputs([1, 2], systems, source, (i, j) -> nothing, i -> [1.0 0], 2)
    @test ynew[1] == [1.0 0]
    @test ynew[2] == [2.0 0]
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
        @test !isstable(multirate_linearization(loop_fast_break, [], [x]; op = merge(opb, Dict(g => 1.02gmax))))
        mlg = multirate_linearization(loop_fast_break, [], [x]; op = merge(opb, Dict(g => 0.98gmax)))
        @test maximum(abs, multipliers(mlg)) < 1
        # Larger classes of perturbations give smaller margins
        @test m.disk_margin_dynamic.α <= m.disk_margin_lti_lower_bound.α + 1e-6
        @test m.disk_margin_lti_lower_bound.γmax <= gmax + 1e-6
    end
end

@testset "Non-nested clocks" begin
    # The periods 2dtf and 3dtf have the common period 6dtf
    ka = ShiftIndex(Clock(2dtf))
    kb = ShiftIndex(Clock(3dtf))
    @variables wa(t) zb(t)
    eqs = [
        D(x) ~ -x + Hold(zb) + d
        wa(ka) ~ 0.5 * wa(ka - 1) + 0.5 * Sample(2dtf)(x)
        zb(kb) ~ 0.6 * zb(kb - 1) - 0.3 * Sample(3dtf)(Hold(wa))
    ]
    @named nonnested = System(eqs, t)
    mln = multirate_linearization(nonnested, [d], [x]; op = Dict(x => 0.0, d => 0.0, wa => 0.0, zb => 0.0))
    @test mln.Tc ≈ 6dtf
    @test mln.event_times ≈ dtf .* [0, 2, 3, 4]
    # One period of the equations on the grid with interval dtf, with the state (x, wa, zb)
    function period(v)
        xv, wv, zv = v
        a = exp(-dtf)
        for m in 0:5
            # `Sample(Hold(wa))` has the value of wa that was held before the tick
            wold = wv
            m % 2 == 0 && (wv = 0.5 * wv + 0.5 * xv)
            m % 3 == 0 && (zv = 0.6 * zv - 0.3 * wold)
            xv = a * xv + (1 - a) * zv
        end
        [xv, wv, zv]
    end
    Mref = reduce(hcat, period.(eachcol(Matrix(1.0I, 3, 3))))
    nonzero(λ) = sort(filter(l -> abs(l) > 1e-8, λ), by = abs)
    @test nonzero(multipliers(mln)) ≈ nonzero(eigvals(Mref))
end

@testset "A signal read within its own partition" begin
    @variables pp(t) qq(t) uu(t)
    eqs = [
        pp(kf) ~ 0.25 * pp(kf - 1) + uu(kf)
        qq(kf) ~ Sample(dtf)(Hold(pp))
    ]
    @named selfread = System(eqs, t)
    ops = Dict(pp => 0.0, qq => 0.0, uu => 0.0)
    mls = multirate_linearization(selfread, [uu], [qq, pp]; op = ops)
    Gs = lifted_ss(mls)
    # `Sample(Hold(pp))` has the value of pp that was held before the tick: qq is pp delayed
    # by one tick
    @test Gs.D ≈ [0; 1]
    z = tf("z", dtf)
    wv = exp10.(range(-1, log10(π / dtf), length = 20))
    @test freqresp(Gs.sys, wv) ≈ freqresp([1 / (z - 0.25); z / (z - 0.25)], wv)
    # The assembly into a single system does not represent the delay
    @test_throws "delays it by one tick" named_ss(selfread, [uu], [qq, pp]; op = ops, hybrid = true)
end

@testset "Event schedule" begin
    using ControlSystemsMTK: event_schedule, common_period
    # Instants that coincide within the tolerance across the boundary of the period are one
    # event, irrespective of the start of the period
    clocks = Dict(1 => (1.0, 1 - 5e-9), 2 => (10.0, 0.0))
    for t0 in (0.0, 0.5, 3.0)
        times, eticks, ticks = event_schedule(clocks, 10.0, t0, 1e-8)
        @test length(times) == 10
        @test count(e -> haskey(e, 2), eticks) == 1
        @test all(e -> haskey(e, 1), eticks)
        @test t0 - 1e-7 <= minimum(times) && maximum(times) < t0 + 10
        @test issorted(times) && issorted(ticks[1])
    end
    @test_throws "not commensurate" common_period([1.0, 0.0010000005])
    @test common_period([0.125, 0.375, 0.25]) ≈ 0.75
end

@testset "Margins of given sensitivity functions" begin
    # A periodic model with one event, one scalar signal and the sensitivity function
    # H*inv(z*I - F)*G + E
    function sensitivity_model(F, G, H, E, Tc)
        port = MultiratePort(:s, [0.0], [1], 1)
        MultirateLinearization([F], [G], [H], [E], [[1]], Tc, 0.0, [0.0], [Symbol("x$i") for i in axes(F, 1)], [port], [port], nothing)
    end
    # S = 1.5 - z^-2: with zero skew the disk margin is 1/‖1 - z^-2‖∞ = 1/2, attained where z^2 = -1
    Tc = 0.125
    m = multirate_margins(sensitivity_model([0.0 0; 1 0], [1.0; 0;;], [0.0 -1], [1.5;;], Tc))
    @test m.disk_margin_dynamic.α ≈ 0.5 rtol = 1e-4
    @test m.disk_margin_dynamic.ω0 ≈ π / (2Tc) rtol = 1e-2
    @test m.disk_margin_lti_lower_bound.α ≈ 0.5 rtol = 1e-4
    # A loss of stability immediately before the gain at which the loop is not well posed
    m = multirate_margins(sensitivity_model([0.0;;], [1e-4;;], [1.0;;], [2.0;;], Tc))
    @test m.gain_interval[2] ≈ (2 + 1e-4) / (1 + 1e-4) rtol = 1e-9
    # A static sensitivity function S = 3 corresponds to the loop-transfer function -2/3,
    # for which the loop is not well posed at the gain 3/2
    m = multirate_margins(sensitivity_model([0.0;;], [0.0;;], [0.0;;], [3.0;;], Tc))
    @test m.gain_interval == (0.0, 1.5)
    # A feedthrough matrix with poor conditioning is well posed for all gains
    port = MultiratePort(:s, [0.0, Tc / 2], [1, 2], 1)
    mlc = MultirateLinearization([zeros(0, 0), zeros(0, 0)], [zeros(0, 2), zeros(0, 2)], [zeros(1, 0), zeros(1, 0)],
        [[1.0 0], [1e7 1]], [[1], [2]], Tc, 0.0, [0.0, Tc / 2], Symbol[], [port], [port], nothing)
    @test multirate_margins(mlc).gain_interval == (0.0, Inf)
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
