using ControlSystemsMTK, ControlSystemsBase, ModelingToolkit, Test
using RobustAndOptimalControl: NamedStateSpace
using ModelingToolkit: t_nounits as t, D_nounits as D, connect

# A stand-in for a discrete compilation pass such as `SynchToolkit.compile_lustre`
discrete_pass(args...) = error("The discrete compilation pass is not expected to be called")
ModelingToolkit.discrete_compile_pass(::typeof(discrete_pass)) = true

function freqresp_isapprox(G1, G2; w = exp10.(range(-2, 1, length = 30)), atol = 1.0e-8)
    isapprox(freqresp(G1, w), freqresp(G2, w); atol)
end

dt = 0.1
k = ShiftIndex(Clock(dt))
@variables x(t) y(t) d(t) r(t) e(t) xi(t) ud(t)
@variables uh(t) [output = true] u(t) [input = true]
kpv, kiv = 2.0, 1.0
@parameters kp = kpv ki = kiv

# A continuous-time first-order plant controlled by a clocked PI controller
controller = [
    e ~ r - Sample(dt)(y)
    xi(k) ~ xi(k - 1) + ki * dt * e(k)
    ud ~ kp * e + xi
]
plant = [
    D(x) ~ -x + u + d
    y ~ x
]
interconnection = [uh ~ Hold(ud); connect(uh, :plant_input, u)]
# The reference and the disturbance are free in `loop`, the reference is fixed in
# `loop_d`, and both are fixed in `loop_closed`
@named loop = System([controller; interconnection; plant], t)
@named loop_d = System([controller; r ~ 0; interconnection; plant], t)
@named loop_closed = System([controller; r ~ 0; d ~ 0; interconnection; plant], t)
op = Dict(x => 0.0, d => 0.0, r => 0.0, xi => 0.0)

z = tf("z", dt)
Pd = c2d(ss(-1.0, 1.0, 1.0, 0.0), dt)
Cd = ss(kpv + kiv * dt * z / (z - 1))
Sd = feedback(1, Pd * Cd)

if isdefined(ModelingToolkit, :linearize_hybrid)
    @testset "named_ss" begin
        G = named_ss(loop_d, [d], [y]; op)
        @test G isa NamedStateSpace
        @test G.Ts == dt
        @test G.u == [Symbol("d(t)")]
        @test G.y == [Symbol("y(t)")]
        @test freqresp_isapprox(G, Pd * Sd)
        @test G.extra[:hybrid_linearization] isa ModelingToolkit.HybridLinearization

        # A discrete compilation pass selects the hybrid linearization
        G = named_ss(loop_d, [d], [y]; op, additional_passes = [discrete_pass])
        @test freqresp_isapprox(G, Pd * Sd)
        @test_throws ModelingToolkit.HybridSystemNotSupportedException named_ss(loop_d, [d], [y]; op, hybrid = false)

        # The result of `ss` is assembled in the same way
        @test freqresp_isapprox(ss(loop_d, [d], [y]; op), Pd * Sd)

        # An analysis point as input
        G = named_ss(loop_closed, loop_closed.plant_input, [y]; op)
        @test G.u == [:loop_closed₊plant_input]
        @test freqresp_isapprox(G, Pd * Sd)

        # With the loop opened, the plant alone remains
        G = named_ss(loop_closed, loop_closed.plant_input, [y]; op, loop_openings = [loop_closed.plant_input])
        @test freqresp_isapprox(G, Pd)
    end

    @testset "Inputs and outputs in several partitions" begin
        # The order of the signals differs from the order of the partitions
        # Without a discrete compilation pass, a model with inputs in a clocked partition is
        # not recognized as hybrid, which is why `hybrid = true` is required
        G = named_ss(loop, [r, d], [ud, y]; op, hybrid = true)
        @test G.u == [Symbol("r(t)"), Symbol("d(t)")]
        @test G.y == [Symbol("ud(t)"), Symbol("y(t)")]
        Gref = [Cd*Sd -Cd*Pd*Sd; Pd*Cd*Sd Pd*Sd]
        @test freqresp_isapprox(G, Gref)
        G = named_ss(loop, [d, r], [y, ud]; op, additional_passes = [discrete_pass])
        @test freqresp_isapprox(G, Gref[[2, 1], [2, 1]])
    end

    @testset "Sensitivity functions" begin
        S = get_named_sensitivity(loop_closed, loop_closed.plant_input; op)
        @test S.Ts == dt
        @test S.u == S.y == [:loop_closed₊plant_input]
        @test freqresp_isapprox(S, Sd)
        T = get_named_comp_sensitivity(loop_closed, loop_closed.plant_input; op)
        @test freqresp_isapprox(T, 1 - Sd)
        L = get_named_looptransfer(loop_closed, loop_closed.plant_input; op)
        @test freqresp_isapprox(L, -Cd * Pd; w = exp10.(range(-1, 1, length = 30)), atol = 1.0e-6)
    end

    @testset "Unsupported models" begin
        dt2 = 0.3
        k2 = ShiftIndex(Clock(dt2))
        @variables w(t) v(t)
        eqs = [
            D(x) ~ -x + Hold(w) + Hold(v) + d
            w(k) ~ 0.5 * w(k - 1) - Sample(dt)(x)
            v(k2) ~ 0.5 * v(k2 - 1) - Sample(dt2)(x)
        ]
        @named multirate = System(eqs, t)
        @test_throws "same periodic clock" named_ss(multirate, [d], [x]; op = Dict(x => 0.0, d => 0.0, w => 0.0, v => 0.0))

        # Direct feedthrough from a held signal to a sampled signal
        eqs = [
            D(x) ~ -x + Hold(ud) + d
            y ~ x + Hold(ud)
            ud ~ -0.5 * Sample(dt)(y)
        ]
        @named feedthrough = System(eqs, t)
        @test_throws "direct feedthrough" named_ss(feedthrough, [d], [x]; op = Dict(x => 0.0, d => 0.0))
    end
else
    @testset "Hybrid linearization unavailable" begin
        @test_throws "linearize_hybrid" named_ss(loop_d, [d], [y]; op)
        @test_throws "linearize_hybrid" named_ss(loop_d, [d], [y]; op, additional_passes = [discrete_pass])
        @test_throws "linearize_hybrid" get_named_sensitivity(loop_closed, loop_closed.plant_input; op)
    end
end

@testset "Functions without support for hybrid models" begin
    @test_throws "clocked partitions" batch_ss(loop_d, [d], [y], [op]; additional_passes = [discrete_pass])
    @test_throws "clocked partitions" build_quadratic_cost_matrix(loop_d, [d], [y => 1.0]; op, additional_passes = [discrete_pass])
end
