# Validation of `multirate_linearization` against simulations of models compiled with
# SynchToolkit: the monodromy matrix is compared with the Jacobian of the map over one common
# period of the clocks, obtained by differentiating the solution of the nonlinear model with
# finite differences.
#
# SynchToolkit is registered in the DyadRegistry only, so this file is not part of the tests
# run in CI. Run it in an environment that provides ControlSystemsMTK, ModelingToolkit,
# SynchToolkit, OrdinaryDiffEqTsit5 and OrdinaryDiffEqNonlinearSolve.
#
# The clock periods are exactly representable in binary floating point. With periods such as
# 0.1 and 0.3, the tick times of the clocks do not coincide exactly in the simulation, and the
# partitions are then evaluated in the order of their tick times rather than simultaneously.
using ModelingToolkit, SynchToolkit, OrdinaryDiffEqTsit5, OrdinaryDiffEqNonlinearSolve
using ControlSystemsMTK, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
A nonlinear plant, a filter on a fast clock and a nonlinear controller on a slow clock that
reads the filter with `Latest`. With `feedthrough != 0`, the sampled signal depends directly
on the held signal. The model has the equilibrium `xc = w = 1`, `z = 0`, at which the
Jacobians of the plant and of the controller differ from those at the origin.

Returns the deviations between the Jacobian of the solution after 1, 2 and 3 periods with
respect to the initial state, and the monodromy matrix raised to the corresponding power.
"""
function monodromy_errors(dtfast, ratio; feedthrough, δ = 1e-4)
    dtslow = ratio * dtfast
    kfast = ShiftIndex(Clock(dtfast))
    kslow = ShiftIndex(Clock(dtslow))
    @variables xc(t) w(t) z(t) y(t)
    @parameters w0 z0
    eqs = [
        D(xc) ~ -xc^3 + Hold(z) + 1
        y ~ xc + feedthrough * Hold(z)
        w(kfast) ~ 0.5w(kfast - 1) + 0.5Sample(dtfast)(y)
        z(kslow) ~ 0.9z(kslow - 1) - 0.2 * (SynchToolkit.Latest()(w, 0.0)^2 - 1)
    ]
    @named sys = System(eqs, t, [xc, w, z, y], [w0, z0];
        initialization_eqs = [w(kfast - 1) ~ w0, z(kslow - 1) ~ z0])
    ssys = mtkcompile(sys; additional_passes = [SynchToolkit.compile_lustre])
    equilibrium = [1.0, 1.0, 0.0, 0.0] # xc, the history of w and z, and the held z
    nperiods = 3
    # The state of the simulated model immediately before the events at the end of each
    # period: xc, the latest w and z, and the held z
    offset = 1e-9
    function states(initial)
        prob = ODEProblem(ssys, [xc => initial[1], w0 => initial[2], z0 => initial[3], Hold(z) => initial[4]],
            (0.0, nperiods * dtslow))
        sol = solve(prob, Tsit5(); abstol = 1e-13, reltol = 1e-13)
        [[sol(m * dtslow - offset; idxs = v) for v in (xc, w, z, Hold(z))] for m in 1:nperiods]
    end
    # The solution remains in the equilibrium
    @test all(s -> s ≈ equilibrium, states(equilibrium))
    # Jacobians of the state after each period with respect to the initial state, by central
    # differences
    perturbed = map(1:4) do j
        e = zeros(4); e[j] = δ
        states(equilibrium + e), states(equilibrium - e)
    end
    J = [reduce(hcat, [(perturbed[j][1][m] - perturbed[j][2][m]) / 2δ for j in 1:4]) for m in 1:nperiods]

    ml = multirate_linearization(sys, [], [xc];
        op = Dict(xc => 1.0, w => 1.0, z => 0.0, w0 => 1.0, z0 => 0.0), warn_missing_op = false)
    F = ControlSystemsMTK.lifted_matrices(ml)[1]
    @test maximum(abs, multipliers(ml)) < 1
    # The state of the periodic model consists of xc, the history variables of w and z, and
    # the held values of z and w. The simulated values correspond to these as follows.
    simulated = map(ml.x) do name
        n = string(name)
        occursin("xc", n) ? 1 : occursin("hold(z", n) ? 4 : occursin("w", n) ? 2 : 3
    end
    # The initial state perturbs xc, the history variables and the held z independently. The
    # held w is not perturbed: it is overwritten at the first event, at which both clocks tick.
    initial = map(ml.x) do name
        n = string(name)
        occursin("hold(z", n) ? 4 : occursin("hold", n) ? 0 : occursin("xc", n) ? 1 : occursin("w", n) ? 2 : 3
    end
    X0 = [Float64(i == j) for i in initial, j in 1:4]
    map(1:nperiods) do m
        norm(J[m][simulated, :] - F^m * X0)
    end
end

@testset "Monodromy matrix and differentiated solution" begin
    for (dtfast, ratio, feedthrough) in ((0.125, 3, 0.0), (0.125, 3, 0.5), (0.25, 2, 0.5), (0.125, 4, -0.3))
        errors = monodromy_errors(dtfast, ratio; feedthrough)
        @info "Deviation of the Jacobians" dtfast ratio feedthrough maximum(errors)
        @test all(<(1e-6), errors)
    end
end
