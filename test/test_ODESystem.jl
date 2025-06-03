using ControlSystemsMTK,
    ControlSystemsBase, ModelingToolkit, RobustAndOptimalControl, Test
import ModelingToolkitStandardLibrary.Blocks as Blocks
using OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqRosenbrock
using LinearAlgebra
conn = ModelingToolkit.connect
connect = ModelingToolkit.connect
## Test SISO (single input, single output) system
@parameters t

P0 = tf(1.0, [1, 1]) |> ss
C0 = pid(1, 1) * tf(1, [0.01, 1]) |> ss

@named P = System(P0)
@test P isa System
# @test length(ModelingToolkit.outputs(P)) == P0.ny
# @test length(ModelingToolkit.inputs(P)) == P0.nu
# @named nonlinear_P = sconnect(x->sign(x)*sqrt(abs(x)), P) # apply input-nonlinearity
@named C        = System(C0)

Pc = complete(P)
op = Dict(Pc.input.u => 0.0)
Q = ControlSystemsMTK.build_quadratic_cost_matrix(Pc, [Pc.input.u], [Pc.x[1] => 2.0]; op)
@test Q[] ≈ 2.0

Q = ControlSystemsMTK.build_quadratic_cost_matrix(Pc, [Pc.input.u], [Pc.output.u => 2.0]; op)
@test Q[] ≈ 2.0

#Mix states and outputs
Q = ControlSystemsMTK.build_quadratic_cost_matrix(Pc, [Pc.input.u], [Pc.x[1] => 2.0, Pc.output.u => 3]; op)
@test Q[] ≈ 2.0 + 3

matrices, ssys = linearize(Pc, [Pc.input.u], [Pc.output.u]; op)

Q = ControlSystemsMTK.build_quadratic_cost_matrix(matrices, ssys, [Pc.x[1] => 2.0])
@test Q[] ≈ 2.0

Q = ControlSystemsMTK.build_quadratic_cost_matrix(matrices, ssys, [Pc.output.u => 2.0])
@test Q[] ≈ 2.0

P1 = ss(Pc, [Pc.input.u], [Pc.output.u]; op)
@test P1 == P0



# === Go the other way, System -> StateSpace ================================
x = unknowns(P) # I haven't figured out a good way to access states, so this is a bit manual and ugly
@unpack input, output = P
P02_named = named_ss(P, [input.u], [output.u]; op)
@test P02_named.x == [Symbol("(x(t))[1]")]
@test P02_named.u == [Symbol("input₊u(t)")]
@test P02_named.y == [Symbol("output₊u(t)")]

P02 = ss(P02_named)
@test P02 == P0 # verify that we get back what we started with

# same for controller
x = unknowns(C)
@nonamespace C02 = named_ss(C, [C.input], [C.output]; op)
@test ss(C02) == C0


## Back again for a complete round trip, test that System get correct names
@named P2 = System(P02_named)
# @test Set(unknowns(P)) ⊆ Set(unknowns(P2))
# @test Set(ModelingToolkit.inputs(P)) ⊆ Set(ModelingToolkit.inputs(P2))
# @test Set(ModelingToolkit.outputs(P)) ⊆ Set(ModelingToolkit.outputs(P2))




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

@named RF = System(tf(1, [0.001, 1])) # ref filter
@named RFv = System(tf(1, [0.001, 1]))
@named P = System(P0) # system model
@named Fv = System(inv(P0i[2, 1])) # feedforward vel
@named Fp = System(inv(P0i[1, 1])) # feedforward pos
@named Cv = System(Cv0) # vel controller
@named Cp = System(Cp0) # pos controller
fb =
    let ref0 = Blocks.Sine(amplitude = 0.2, frequency = 1, name = :r),
        disturbance = Blocks.Step(height = 100, start_time = 10, name = :d)

        @named input = Blocks.RealInput()
        @named output = Blocks.RealOutput()
        Dₜ = Differential(t)
        fb = System(
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
simplified_sys = mtkcompile(fb)


x0 = Pair[
    P.x[1] => 0.0
    P.x[3] => 0.0
]

prob = ODEProblem(simplified_sys, x0, (0.0, 20.0), jac = true)
sol = solve(prob, Rodas5(), saveat = 0:0.01:20)
if isinteractive()
    @show sol.retcode
    plot(sol, layout = length(unknowns(simplified_sys)) + 1)
    plot!(sol.t, sol[P.x[1]] - sol[P.x[3]], sp = 12, lab = "Δq")

    ##
    plot(sol.t, sol[P.x[1]] - sol[fb.r.output.u], lab = "qₘ", title = "Control error")
    plot!(sol.t, sol[P.x[3]] - sol[fb.r.output.u], lab = "qₐ")
end

## Double-mass model in MTK
using ModelingToolkit, OrdinaryDiffEqRosenbrock, LinearAlgebra
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

@named inertia1 = Inertia(; J = m1, w=0)
@named inertia2 = Inertia(; J = m2, w=0)

@named spring = Spring(; c = k)
@named damper = Damper(; d = c)

@named torque = Torque(use_support = false)

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
        return @named model = System(eqs, t; systems = [sens, torque, inertia1, inertia2, spring, damper, u])
    end
    System(eqs, t; systems = [sens, torque, inertia1, inertia2, spring, damper], name)
end

model = SystemModel() |> complete

op = Dict(model.inertia1.flange_b.phi => 0.0, model.torque.tau.u => 0)
lsys = named_ss(model, [model.torque.tau.u], [model.inertia1.phi, model.inertia2.phi]; op)
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
defs = merge(Dict(unknowns(model) .=> 0), defs)
_, p = ModelingToolkit.get_u0_p(ssys, defs, defs)


sympars = ModelingToolkit.parameters(ssys)

fun = Symbolics.build_function(sys, sympars; expression=Val{false}, force_SA=true)
fun(p)
@test @allocated(fun(p)) <= 256
static_lsys = fun(p)

@test static_lsys == lsys.sys


## Named sensitivity funcitons

@named P = Blocks.FirstOrder(k = 1, T = 1)
@named C = Blocks.Gain(; k = 1)
@named add = Blocks.Add(k2 = -1)
t = ModelingToolkit.get_iv(P)

eqs = [connect(P.output, :plant_output, add.input2)
    connect(add.output, C.input)
    connect(C.output, :plant_input, P.input)]

sys_inner = System(eqs, t, systems = [P, C, add], name = :inner)

@named r = Blocks.Constant(k = 1)
@named F = Blocks.FirstOrder(k = 1, T = 3)

eqs = [connect(r.output, F.input)
    connect(F.output, sys_inner.add.input1)]
sys_outer = System(eqs, t, systems = [F, sys_inner, r], name = :outer)

matrices, _ = Blocks.get_sensitivity(sys_outer, [sys_outer.inner.plant_input, sys_outer.inner.plant_output])
S = ss(matrices...)

Sn = get_named_sensitivity(sys_outer, [sys_outer.inner.plant_input, sys_outer.inner.plant_output])

@test S == Sn.sys

@test Sn.u == Sn.y == [:outer₊inner₊plant_input, :outer₊inner₊plant_output] == [:outer₊inner₊plant_input, :outer₊inner₊plant_output]


## Test connector names
P = named_ss(ssrand(1,1,1), u=:jörgen, y=:solis)
@named Pode = System(P)
ModelingToolkit.isconnector(Pode.jörgen)
ModelingToolkit.isconnector(Pode.solis)



## Test causal simplification
using LinearAlgebra
using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition
using Test

using ControlSystemsMTK
using ControlSystemsMTK.ControlSystemsBase: sminreal, minreal, poles, ss, tf
connect = ModelingToolkit.connect

@independent_variables t
D = Differential(t)

@named link1 = Link(; m = 0.2, l = 10, I = 1, g = -9.807)
@named cart = TranslationalPosition.Mass(; m = 1, s = 0)
@named fixed = TranslationalPosition.Fixed()
@named force = Force(use_support = false)

eqs = [connect(link1.TX1, cart.flange)
       connect(cart.flange, force.flange)
       connect(link1.TY1, fixed.flange)]

@named model = System(eqs, t, [], []; systems = [link1, cart, force, fixed])
lin_outputs = [cart.s, cart.v, link1.A, link1.dA]
lin_inputs = [force.f.u]

# => nothing to remove extra defaults
op = Dict(cart.s => 10, cart.v => 0, link1.A => -pi/2, link1.dA => 0, force.f.u => 0, link1.x1 => nothing, link1.y1 => nothing, link1.x2 => nothing, link1.x_cm => nothing)
G = named_ss(model, lin_inputs, lin_outputs; allow_symbolic = true, op,
    allow_input_derivatives = true, zero_dummy_der = true)
G = sminreal(G)
@info "minreal"
G = minreal(G)
@info "poles"
ps = poles(G)

@test minimum(abs, ps) < 1e-6
@test minimum(abs, complex(0, 1.3777260367206716) .- ps) < 1e-10

lsys, syss = linearize(model, lin_inputs, lin_outputs, op = op,
    allow_input_derivatives = true)
lsyss, sysss = ModelingToolkit.linearize_symbolic(model, lin_inputs, lin_outputs;
    allow_input_derivatives = true)

dummyder = setdiff(unknowns(sysss), unknowns(model))
# op2 = merge(ModelingToolkit.guesses(model), op, Dict(x => 0.0 for x in dummyder))
op2 = merge(ModelingToolkit.defaults(syss), op)
op2[link1.fy1] = -op2[link1.g] * op2[link1.m]
op2[cart.f] = 0

@test substitute(lsyss.A, op2) ≈ lsys.A
# We cannot pivot symbolically, so the part where a linear solve is required
# is not reliable.
@test substitute(lsyss.B, op2)[1:6, 1] ≈ lsys.B[1:6, 1]
@test substitute(lsyss.C, op2) ≈ lsys.C
@test substitute(lsyss.D, op2) ≈ lsys.D

@test G.nx == 4
@test G.nu == length(lin_inputs)
@test G.ny == length(lin_outputs)

## Test difficult `named_ss` simplification
using ControlSystemsMTK, ControlSystemsBase, RobustAndOptimalControl
lsys = (A = [0.0 2.778983834717109e8 1.4122312296634873e6 0.0; 0.0 0.0 0.0 0.037848975765016724; 0.0 24.837541148074962 0.12622006230897712 0.0; -0.0 -4.620724819774693 -0.023481719514324866 -0.6841991610512456], B = [-5.042589978197361e8 0.0; -0.0 0.0; -45.068824982602656 -0.0; 8.384511049369085 54.98555939873381], C = [0.0 0.0 0.954929658551372 0.0], D = [0.0 0.0])

# lsys = (A = [-0.0075449237853825925 1.6716817118020731e-6 0.0; 1864.7356343162514 -0.4131578457122937 0.0; 0.011864343330426718 -2.6287085638214332e-6 0.0], B = [0.0 0.0; 0.0 52566.418015009294; 0.0 0.3284546792274811], C = [1.4683007399899215e8 0.0 0.0], D = [-9.157636303058283e7 0.0])

G = ControlSystemsMTK.causal_simplification(lsys, [1=>2])
G2 = ControlSystemsMTK.causal_simplification(lsys, [1=>2], descriptor=false)
G2 = minreal(G2, 1e-12)

@test dcgain(G, 1e-5)[] ≈ dcgain(G2, 1e-5)[] rtol=1e-3
@test freqresp(G, 1)[] ≈ freqresp(G2, 1)[]
@test freqresp(G, 10)[] ≈ freqresp(G2, 10)[]

z = 0.462726166562343204837317130554462562

@test minimum(abs, tzeros(G) .- z) < sqrt(eps())
@test minimum(abs, tzeros(G2) .- z) < sqrt(eps())

# using GenericSchur

Gb = balance_statespace(G)[1]
Gb = minreal(Gb, 1e-8)
@test Gb.nx == 2
@test minimum(abs, tzeros(Gb) .- z) < sqrt(eps())

w = exp10.(LinRange(-12, 2, 2000))
# ControlSystemsBase.bodeplot([G, G2, minreal(G, 1e-8)], w)


##

# S = schur(A,B)
# V = eigvecs(S)
