# Validation of `multirate_linearization` against simulations of models compiled with
# SynchToolkit. SynchToolkit is registered in the DyadRegistry only, so this file is not part
# of the tests run in CI. Run it in an environment that provides ControlSystemsMTK,
# ModelingToolkit, SynchToolkit, OrdinaryDiffEqTsit5 and OrdinaryDiffEqNonlinearSolve.
#
# The clock periods are exactly representable in binary floating point. With periods such as
# 0.1 and 0.3, the tick times of the clocks do not coincide exactly in the simulation, and the
# partitions are then evaluated in the order of their tick times rather than simultaneously.
using ModelingToolkit, SynchToolkit, OrdinaryDiffEqTsit5, OrdinaryDiffEqNonlinearSolve
using ControlSystemsMTK, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

# The largest deviation, over a number of periods, between the state of the simulated model
# and the state propagated with the monodromy matrix
function monodromy_error(dtfast, ratio; feedthrough)
    dtslow = ratio * dtfast
    kfast = ShiftIndex(Clock(dtfast))
    kslow = ShiftIndex(Clock(dtslow))
    @variables xc(t) w(t) z(t) y(t)
    eqs = [
        D(xc) ~ -xc + Hold(z)
        # With `feedthrough != 0`, the sampled signal depends directly on the held signal
        y ~ xc + feedthrough * Hold(z)
        w(kfast) ~ 0.5w(kfast - 1) + 0.5Sample(dtfast)(y)
        # The slow partition reads the fast partition at coinciding ticks
        z(kslow) ~ 0.9z(kslow - 1) - 0.4 * SynchToolkit.Latest()(w, 0.0)
    ]
    @named sys = System(eqs, t; initialization_eqs = [w(kfast - 1) ~ 0.2, z(kslow - 1) ~ -0.3])
    ssys = mtkcompile(sys; additional_passes = [SynchToolkit.compile_lustre])
    nperiods = 5
    prob = ODEProblem(ssys, [xc => 1.0], (0.0, nperiods * dtslow + dtslow / 3))
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
    ml = multirate_linearization(sys, [], [xc]; op = Dict(xc => 0.0, w => 0.0, z => 0.0), warn_missing_op = false)
    F = ControlSystemsMTK.lifted_matrices(ml)[1]
    # The state before the events at `tt`: `xc` and the latest values of `w` and `z`, which
    # are both the history variables and the held values
    offset = 1e-9
    state(tt) = map(ml.x) do name
        n = string(name)
        v = occursin("xc", n) ? xc : occursin("w", n) ? w : z
        sol(tt - offset; idxs = v)
    end
    xk = state(dtslow)
    maximum(2:nperiods) do m
        xk = F * xk
        norm(xk - state(m * dtslow))
    end
end

@testset "Monodromy matrix and simulation" begin
    @test monodromy_error(0.125, 3; feedthrough = 0.0) < 1e-7
    @test monodromy_error(0.125, 3; feedthrough = 0.5) < 1e-7
    @test monodromy_error(0.25, 2; feedthrough = 0.5) < 1e-7
end
