using ControlSystemsMTK,
    ControlSystemsBase, ModelingToolkit, OrdinaryDiffEq, RobustAndOptimalControl
import ModelingToolkitStandardLibrary.Blocks as Blocks
conn = ModelingToolkit.connect
## Test SISO (single input, single output) system
@parameters t

P0 = tf(1.0, [1, 1]) |> ss
C0 = pid(1, 1) * tf(1, [0.01, 1]) |> ss

@named P = ODESystem(P0)
@test P isa ODESystem
@test length(ModelingToolkit.outputs(P)) == P0.ny
@test length(ModelingToolkit.inputs(P)) == P0.nu
# @named nonlinear_P = sconnect(x->sign(x)*sqrt(abs(x)), P) # apply input-nonlinearity
@named C        = ODESystem(C0)
@named loopgain = sconnect(C, P)
@named ref      = Blocks.Sine(frequency = 1)
fb              = feedback(loopgain, name = :fb) * ref
fb              = structural_simplify(fb)

@test length(states(P)) == 3 # 1 + u + y
@test length(states(C)) == 4 # 2 + u + y

x0 = Pair[loopgain.P.x[1]=>1]

prob = ODEProblem(fb, x0, (0.0, 10.0))
sol = solve(prob, Rodas5())
@test Int(sol.retcode) == 1
isinteractive() && plot(sol)

fb              = structural_simplify(sconnect(sin, P))
prob = ODEProblem(fb, x0, (0.0, 10.0))
sol2 = solve(prob, Rodas5())
@test Int(sol2.retcode) == 1
isinteractive() && plot(sol2)

Pc = complete(P)
Q = ControlSystemsMTK.build_quadratic_cost_matrix(Pc, [Pc.input.u], [Pc.x[1] => 2.0])
@test Q[] ≈ 2.0

Q = ControlSystemsMTK.build_quadratic_cost_matrix(Pc, [Pc.input.u], [Pc.output.u => 2.0])
@test Q[] ≈ 2.0

#Mix states and outputs
Q = ControlSystemsMTK.build_quadratic_cost_matrix(Pc, [Pc.input.u], [Pc.x[1] => 2.0, Pc.output.u => 3])
@test Q[] ≈ 2.0 + 3

matrices, ssys = linearize(Pc, [Pc.input.u], [Pc.output.u])

Q = ControlSystemsMTK.build_quadratic_cost_matrix(matrices, ssys, [Pc.x[1] => 2.0])
@test Q[] ≈ 2.0

Q = ControlSystemsMTK.build_quadratic_cost_matrix(matrices, ssys, [Pc.output.u => 2.0])
@test Q[] ≈ 2.0

P1 = ss(Pc, [Pc.input.u], [Pc.output.u])
@test P1 == P0



# === Go the other way, ODESystem -> StateSpace ================================
x = states(P) # I haven't figured out a good way to access states, so this is a bit manual and ugly
@unpack input, output = P
P02_named = named_ss(P, [input.u], [output.u])
@test P02_named.x == [Symbol("(x(t))[1]")]
@test P02_named.u == [Symbol("input₊u(t)")]
@test P02_named.y == [Symbol("output₊u(t)")]

P02 = ss(P02_named)
@test P02 == P0 # verify that we get back what we started with

# same for controller
x = states(C)
@nonamespace C02 = named_ss(C, [C.input], [C.output])
@test ss(C02) == C0


## Back again for a complete round trip, test that ODESystem get correct names
@named P2 = ODESystem(P02_named)
@test Set(states(P2)) == Set(states(P))
@test Set(ModelingToolkit.inputs(P2)) == Set(ModelingToolkit.inputs(P))
@test Set(ModelingToolkit.outputs(P2)) == Set(ModelingToolkit.outputs(P))




## Test more complicated feedback connection =====================================================
@parameters t

Jm = 10
Ja = 30
k  = 4e4
c0 = 100
c1 = 100
c2 = 0.001
tx = 20

A = [
    0.0 1 0 0
    -k/Jm -(c1 + c0)/Jm k/Jm c1/Jm
    0 0 0 1
    k/Ja c1/Ja -k/Ja -(c1 + c2)/Ja
]

B  = [0, 1 / Jm, 0, 0]
Bd = [0, 0, 0, 1 / Ja]
B  = [B Bd]
C  = [1 0 0 0
      0 1 0 0]

isinteractive() && @info "Resonance frequency: $((imag.(eigvals(A)) ./ (2π)))"

#

P0 = ss(A, B, C, 0)

s = tf("s")
Cv0 = tx * pid(4, 2) * tf(1, [0.01, 1])  #|> ss
Cp0 = tx * pid(2, 0) * tf(1, [0.01, 1])  #|> ss

P0i = deepcopy(P0)
P0i.D .= 1e-8

@named RF = ODESystem(tf(1, [0.001, 1])) # ref filter
@named RFv = ODESystem(tf(1, [0.001, 1]))
@named P = ODESystem(P0) # system model
@named Fv = ODESystem(inv(P0i[2, 1])) # feedforward vel
@named Fp = ODESystem(inv(P0i[1, 1])) # feedforward pos
@named Cv = ODESystem(Cv0) # vel controller
@named Cp = ODESystem(Cp0) # pos controller
fb =
    let ref0 = Blocks.Sine(amplitude = 0.2, frequency = 1, name = :r),
        disturbance = Blocks.Step(height = 100, start_time = 10, name = :d)

        @named input = Blocks.RealInput()
        @named output = Blocks.RealOutput()
        Dₜ = Differential(t)
        fb = ODESystem(
            [
                conn(ref0.output, RF.input) # filter position reference
                expand_derivatives(Dₜ(ref0.output.u) ~ RFv.input.u) # Filter vel reference
                RF.output.u - P.output.u[1] ~ Cp.input.u # pos controller input is pos error
                disturbance.output.u ~ P.input.u[2] # disturbance enters on arm acceleration
                Cp.output.u + RFv.output.u - P.output.u[2] ~ Cv.input.u # vel controller input is vel error + pos controller output
                Cv.output.u + Fp.output.u ~ P.input.u[1] # robot input is vel controller output and trq ff
                output.u ~ P.output.u[1] # output is robot motor pos
                RF.output.u ~ Fp.input.u
            ],
            t;
            systems = [P, Cv, Cp, Fp, RF, RFv, input, output, ref0, disturbance],
            name = :feedback,
        )
    end
simplified_sys = structural_simplify(fb)


x0 = Pair[
    P.x[1] => 0.0
    P.x[3] => 0.0
]
p = Pair[]

prob = ODEProblem(simplified_sys, x0, (0.0, 20.0), p, jac = true)
sol = solve(prob, OrdinaryDiffEq.Rodas5(), saveat = 0:0.01:20)
if isinteractive()
    @show sol.retcode
    plot(sol, layout = length(states(simplified_sys)) + 1)
    plot!(sol.t, sol[P.x[1]] - sol[P.x[3]], sp = 12, lab = "Δq")

    ##
    plot(sol.t, sol[P.x[1]] - sol[fb.r.output.u], lab = "qₘ", title = "Control error")
    plot!(sol.t, sol[P.x[3]] - sol[fb.r.output.u], lab = "qₐ")
end

## Double-mass model in MTK
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks: Sine
using ModelingToolkit: connect
import ModelingToolkitStandardLibrary.Blocks
t = Blocks.t

# Parameters
m1 = 1
m2 = 1
k = 1000 # Spring stiffness
c = 10   # Damping coefficient

@named inertia1 = Inertia(; J = m1)
@named inertia2 = Inertia(; J = m2)

@named spring = Spring(; c = k)
@named damper = Damper(; d = c)

@named torque = Torque()

function SystemModel(u=nothing; name=:model)
    @named sens = Rotational.AngleSensor()
    eqs = [
        connect(torque.flange, inertia1.flange_a)
        connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
        connect(inertia2.flange_a, spring.flange_b, damper.flange_b)
        connect(inertia2.flange_b, sens.flange)
    ]
    if u !== nothing 
        push!(eqs, connect(u.output, :u, torque.tau))
        return @named model = ODESystem(eqs, t; systems = [sens, torque, inertia1, inertia2, spring, damper, u])
    end
    ODESystem(eqs, t; systems = [sens, torque, inertia1, inertia2, spring, damper], name)
end

model = SystemModel() |> complete

lsys = named_ss(model, [model.torque.tau.u], [model.inertia1.phi, model.inertia2.phi])
@test -1000 ∈ lsys.A
@test -10 ∈ lsys.A
@test 1000 ∈ lsys.A
@test 10 ∈ lsys.A
@test 1 ∈ lsys.B
@test 1 ∈ lsys.C

# model = SystemModel(Sine(frequency=30/2pi, name=:u)) |> complete
# lsys = named_ss(model, :u, [model.inertia1.phi, model.inertia2.phi])

##
mats, ssys = ModelingToolkit.linearize_symbolic(model, [model.torque.tau.u], [model.inertia1.phi, model.inertia2.phi])
sys = ss((mats...,)[1:4]...)


defs = ModelingToolkit.defaults(ssys)
sympars = ModelingToolkit.parameters(ssys)
_, p = ModelingToolkit.get_u0_p(ssys, defs, defs)

fun = Symbolics.build_function(sys, sympars; expression=Val{false}, force_SA=true)
fun(p)
@test @allocated(fun(p)) <= 256
static_lsys = fun(p)

@test static_lsys == lsys.sys
