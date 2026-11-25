# Validation of `multirate_linearization` against simulations of models compiled with
# SynchToolkit: the matrices of the lifted system are compared with Jacobians of the solution
# of the nonlinear model over common periods of the clocks, obtained by finite differences
# with respect to the initial state and to the inputs.
#
# SynchToolkit is registered in the DyadRegistry only, so this file is not part of the tests
# run in CI. Run it in an environment that provides ControlSystemsMTK, ModelingToolkit,
# SynchToolkit, OrdinaryDiffEqTsit5 and OrdinaryDiffEqNonlinearSolve.
#
# The clock periods and phases are exactly representable in binary floating point. With
# periods such as 0.1 and 0.3, the tick times of the clocks do not coincide exactly in the
# simulation, and the partitions are then evaluated in the order of their tick times rather
# than simultaneously.
using ModelingToolkit, SynchToolkit, OrdinaryDiffEqTsit5, OrdinaryDiffEqNonlinearSolve
using ControlSystemsMTK, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
A nonlinear plant with the disturbance input `d`, a filter on a fast clock, and a nonlinear
controller with the reference input `r` on a slow clock with phase `phase` that reads the
filter with `Latest`. The sampled signal `y` depends directly on the held signal with the
gain `feedthrough` and on the disturbance with the gain `input_feedthrough`. The model has
the equilibrium `xc = w = 1`, `z = 0` for `d = 0`, `r = 1`, at which the Jacobians of the
plant and of the controller differ from those at the origin.

The lifted system with the inputs `d`, `r` and the outputs `y`, `w`, `z`, for a period that
starts at `t0`, is compared with the Jacobians of the simulated model with respect to the
state at the start of the period and to the inputs during the period: the outputs during
that period and the next one, and the state after one and two periods. Returns the
deviations.
"""
function lifted_system_errors(dtfast, ratio; feedthrough, input_feedthrough = 0.3, phase = 0.0, t0 = 0.0, δ = 1e-4)
    dtslow = ratio * dtfast
    Tc = dtslow
    fastclock = Clock(dtfast)
    slowclock = Clock(dtslow; phase)
    kfast = ShiftIndex(fastclock)
    kslow = ShiftIndex(slowclock)
    @variables xc(t) w(t) z(t) y(t) d(t) r(t) rc(t)
    @parameters w0 z0
    eqs = [
        D(xc) ~ -xc^3 + Hold(z) + 1 + d
        y ~ xc + feedthrough * Hold(z) + input_feedthrough * d
        w(kfast) ~ 0.5w(kfast - 1) + 0.5Sample(fastclock)(y)
        z(kslow) ~ 0.9z(kslow - 1) - 0.2 * (SynchToolkit.Latest()(w, 0.0)^2 - r)
    ]

    # The linearization, with `d` and `r` as inputs
    @named lin = System(eqs, t)
    ml = multirate_linearization(lin, [d, r], [y, w, z]; t0,
        op = Dict(xc => 1.0, w => 1.0, z => 0.0, d => 0.0, r => 1.0), warn_missing_op = false)
    F, G, H, E = ControlSystemsMTK.lifted_matrices(ml)
    @test maximum(abs, multipliers(ml)) < 1
    dport, rport = ml.inputs
    nd = length(dport.times)
    NU = nd + 1
    # The instants of the signals, derived from the clocks independently of the model: the
    # period starts at a tick of the fast clock, at which d, y and w have their first instant
    slowtick = phase + dtslow * ceil((t0 - phase) / dtslow)
    @test dport.times ≈ t0 .+ dtfast .* (0:ratio-1)
    @test rport.times ≈ [slowtick]
    @test ml.outputs[1].times ≈ dport.times
    @test ml.outputs[2].times ≈ dport.times
    @test ml.outputs[3].times ≈ [slowtick]

    # The simulated model: `d` takes the values `dv` on the intervals that start at the
    # instants of its port, and the reference takes the value `rv` during the period
    @parameters dv[1:nd] rv
    breaks = [dport.times; dport.times[end] + dtfast]
    dexpr = 0.0
    for i in nd:-1:1
        dexpr = ifelse(t < breaks[i+1], dv[i], dexpr)
    end
    dexpr = ifelse(t < breaks[1], 0.0, dexpr)
    simeqs = [
        eqs
        d ~ dexpr
        rc ~ ifelse((t >= t0) & (t < t0 + Tc), rv, 1.0)
        r ~ Sample(slowclock)(rc)
    ]
    @named sim = System(simeqs, t, [xc, w, z, y, d, r, rc], [w0, z0, dv, rv];
        initialization_eqs = [w(kfast - 1) ~ w0, z(kslow - 1) ~ z0])
    ssim = mtkcompile(sim; additional_passes = [SynchToolkit.compile_lustre])
    @test isequal(unknowns(ssim), [xc])

    # Parameters of the experiment: the initial xc, the initial history of w and z, the
    # initially held z, the values of d, and the value of r
    baseline = [1.0; 1.0; 0.0; 0.0; zeros(nd); 1.0]
    offset = 1e-9
    statevars = (xc, w, z, Hold(z))
    outvars = (y, w, z)
    function observe(θ)
        prob = ODEProblem(ssim, [xc => θ[1], w0 => θ[2], z0 => θ[3], Hold(z) => θ[4], dv => θ[5:4+nd], rv => θ[end]],
            (0.0, t0 + 2Tc + dtfast))
        sol = solve(prob, Tsit5(); abstol = 1e-13, reltol = 1e-13, tstops = [breaks; t0; t0 + Tc])
        # The state immediately before the events at `tt`
        state(tt) = [sol(tt - offset; idxs = v) for v in statevars]
        # The outputs immediately after the events at their instants, in the order of the
        # lifted channels
        function outputs(shift)
            yv = zeros(size(H, 1))
            for (port, v) in zip(ml.outputs, outvars), (i, ti) in enumerate(port.times)
                yv[port.channels[i]] = sol(ti + shift + offset; idxs = v)
            end
            yv
        end
        (; s0 = t0 > 0 ? state(t0) : zeros(4), y1 = outputs(0.0), s1 = state(t0 + Tc), y2 = outputs(Tc), s2 = state(t0 + 2Tc))
    end
    # The solution remains in the equilibrium
    nominal = observe(baseline)
    @test nominal.s1 ≈ [1, 1, 0, 0] atol = 1e-10
    @test nominal.s2 ≈ [1, 1, 0, 0] atol = 1e-10
    # Jacobians by central differences
    nθ = length(baseline)
    plus = [observe(baseline + δ * (1:nθ .== j)) for j in 1:nθ]
    minus = [observe(baseline - δ * (1:nθ .== j)) for j in 1:nθ]
    jac(field) = reduce(hcat, [(getfield(plus[j], field) - getfield(minus[j], field)) / 2δ for j in 1:nθ])

    # The state of the periodic model consists of xc, the history variables of w and z, the
    # held values of z and w, and the held value of d. `simulated[i]` is the index among the
    # simulated variables (xc, w, z, held z) that the state variable i equals before an
    # event, and zero for the held d, which is not a variable of the simulated model.
    @test length(ml.x) == 6
    simulated = map(ml.x) do name
        n = string(name)
        n == "hold(d(t))" ? 0 : n == "hold(z(t))" ? 4 : n == "hold(w(t))" ? 2 : n == "xc(t)" ? 1 :
        n == "wₜ₋₁(t)" ? 2 : n == "zₜ₋₁(t)" ? 3 : error("Unexpected state variable $n")
    end
    known = findall(!iszero, simulated)
    # The held d is not a variable of the simulated model. After a period it equals the last
    # value of d in the period, irrespective of the state at the start of the period.
    held_d = only(findall(iszero, simulated))
    @test iszero(F[held_d, :])
    @test G[held_d, :] == (1:NU .== dport.channels[end])
    M = [Float64(simulated[i] == j) for i in eachindex(simulated), j in 1:4]
    # Jacobian of the state at the start of the period
    if t0 > 0
        X0 = M * jac(:s0)
    else
        # The initial values perturb xc, the history variables and the held z. The held w is
        # not an initial value of the simulated model; it is overwritten by the fast
        # partition before it is read.
        initial = map(ml.x) do name
            n = string(name)
            n == "hold(z(t))" ? 4 : startswith(n, "hold") ? 0 : n == "xc(t)" ? 1 : n == "wₜ₋₁(t)" ? 2 : 3
        end
        X0 = [Float64(i == j) for i in initial, j in 1:nθ]
    end
    # The inputs of the period, in the order of the lifted channels
    U = zeros(NU, nθ)
    for i in 1:nd
        U[dport.channels[i], 4+i] = 1
    end
    U[rport.channels[1], nθ] = 1

    X1 = F * X0 + G * U
    (;
        output_first_period = norm(jac(:y1) - (H * X0 + E * U)),
        state_first_period = norm((M*jac(:s1))[known, :] - X1[known, :]),
        output_second_period = norm(jac(:y2) - H * X1),
        state_second_period = norm((M*jac(:s2))[known, :] - (F*X1)[known, :]),
        # The magnitudes of the compared quantities
        scale = (; E = norm(E), G = norm(G), H = norm(H), F = norm(F)),
    )
end

@testset "Lifted system and differentiated solution" begin
    dtfast = 0.125
    cases = [
        (; ratio = 3, feedthrough = 0.0, phase = 0.0, t0 = 0.0),
        (; ratio = 3, feedthrough = 0.5, phase = 0.0, t0 = 0.0),
        (; ratio = 2, feedthrough = 0.5, phase = 0.0, t0 = 0.0),
        (; ratio = 4, feedthrough = -0.3, phase = 0.0, t0 = 0.0),
        # A clock with a phase. The solution of the simulated model cannot be evaluated before
        # the first tick of every clock, so the period starts one common period later.
        # The slow clock ticks together with the second tick of the fast clock
        (; ratio = 3, feedthrough = 0.5, phase = dtfast, t0 = 3dtfast),
        # The slow clock ticks between the ticks of the fast clock
        (; ratio = 3, feedthrough = 0.5, phase = dtfast / 2, t0 = 3dtfast),
        # The period starts at another tick of the fast clock than the one at which the slow
        # clock ticks
        (; ratio = 3, feedthrough = 0.5, phase = 0.0, t0 = 4dtfast),
        (; ratio = 3, feedthrough = 0.5, phase = dtfast / 2, t0 = 5dtfast),
    ]
    for c in cases
        errors = lifted_system_errors(dtfast, c.ratio; c.feedthrough, c.phase, c.t0)
        @info "Deviation of the Jacobians" c errors
        @test errors.output_first_period < 1e-6
        @test errors.state_first_period < 1e-6
        @test errors.output_second_period < 1e-6
        @test errors.state_second_period < 1e-6
    end
end
