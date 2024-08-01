using ControlSystemsMTK, ModelingToolkit, RobustAndOptimalControl, MonteCarloMeasurements
using ModelingToolkit: getdefault
unsafe_comparisons(true)

# Create a model
@parameters t k=10 k3=2 c=1
@variables x(t)=0 [bounds = (-2, 2)]
@variables v(t)=0

D = Differential(t)
@named y = Blocks.RealOutput()
@named u = Blocks.RealInput()

eqs = [D(x) ~ v
       D(v) ~ -k * x - k3 * x^3 - c * v + 10u.u
       y.u ~ x]


@named duffing = ODESystem(eqs, t, systems=[y, u])

bounds = getbounds(duffing, unknowns(duffing))
sample_within_bounds((l, u)) = (u - l) * rand() + l
# Create a vector of operating points
N = 10
xs = range(getbounds(x)[1], getbounds(x)[2], length=N)
ops = Dict.(x .=> xs)

inputs, outputs = [u.u], [y.u]
Ps, ssys = batch_ss(duffing, inputs, outputs , ops)
@test length(Ps) == N

@test Ps[1] == ss(linearize(duffing, inputs, outputs; op=ops[1])[1]...)
@test Ps[end] == ss(linearize(duffing, inputs, outputs; op=ops[end])[1]...)

##

using DataInterpolations
@named Cgs = GainScheduledStateSpace(Ps, xs, interpolator=LinearInterpolation)
@test Cgs isa ODESystem
# This is tested better in the docs

## C-code generation
# using SymbolicControlSystems
# code = SymbolicControlSystems.print_c_array(stdout, Ps, xs, "gain_scheduled_controller", struct_name="hej", struct_type="kaj")

## Simulate
using OrdinaryDiffEq, ControlSystemsBase
import ModelingToolkitStandardLibrary.Blocks

@named fb = Blocks.Add(k2=-1)
@named ref = Blocks.Square(frequency=1/6, amplitude=0.5, offset=0.5, start_time=1)
@named F = Blocks.SecondOrder(w=10, d=0.7)
@named C = ODESystem(pid(1,1,0; state_space=true, Tf=0.01))


closed_loop_eqs = [
    connect(ref.output, :r, F.input)
    connect(F.output, fb.input1)
    connect(duffing.y, :y, fb.input2)
    connect(fb.output, C.input)
    connect(C.output, duffing.u)
]

@named closed_loop = ODESystem(closed_loop_eqs, t, systems=[duffing, C, fb, ref, F])

ssys = structural_simplify(closed_loop)
prob = ODEProblem(ssys, [F.xd => 0], (0.0, 8.0))
sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
# plot(sol)

## Linearize around trajectory

time = 0:0.1:8
inputs, outputs = [duffing.u.u], [duffing.y.u]
Ps2, ssys = trajectory_ss(closed_loop, :r, :y, sol; t=time)
@test length(Ps2) == length(time)
# bodeplot(Ps2)
