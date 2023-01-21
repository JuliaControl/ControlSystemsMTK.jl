# Batch Linearization and gain scheduling
This example will demonstrate how to linearize a nonlinear ModelingToolkit model in multiple different operating points, and some tools to work with groups of linear models representing the same system in different operating points. We'll end with designing and simulating a [gain-scheduled controller](https://en.wikipedia.org/wiki/Gain_scheduling), i.e., a nonlinear controller created as an interpolation between linear controllers.


## System model
The model will be a simple Duffing oscillator:
```@example BATCHLIN
using ControlSystemsMTK, ModelingToolkit, MonteCarloMeasurements, ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: getdefault
unsafe_comparisons(true)

# Create a model
@parameters t k=10 k3=2 c=1
@variables x(t)=0 [bounds = (-0.5, 1.5)]
@variables v(t)=0

D = Differential(t)

@named y = Blocks.RealOutput()
@named u = Blocks.RealInput()

eqs = [D(x) ~ v
       D(v) ~ -k * x - k3 * x^3 - c * v + 10u.u
       y.u ~ x]


@named duffing = ODESystem(eqs, t, systems=[y, u])
```

## Batch linearization
To perform batch linearization, we create a vector of operating points, and then linearize the model around each of these points. The function [`batch_ss`](@ref) does this for us, and returns a vector of `StateSpace` models, one for each operating point. An operating point is a `Dict` that maps variables in the MTK model to numerical values. In the example below, we simply sample the variables uniformly within their bounds specified when we created the variables (normally, we might want to linearize on stationary points)
```@example BATCHLIN
N = 16 # Number of samples
xs = range(getbounds(x)[1], getbounds(x)[2], length=N)
ops = Dict.(x .=> xs)
```

Just like [`ModelingToolkit.linearize`](@ref), [`batch_ss`](@ref) takes the set of inputs and the set of outputs to linearize between.
```@example BATCHLIN
Ps, ssys = batch_ss(duffing, [u.u], [y.u], ops)
nothing # hide
```

Plotting functions like [`bodeplot`](@ref) accept vectors of systems, so this works
```@example BATCHLIN
using ControlSystemsBase, Plots
w = exp10.(LinRange(-2, 3, 200))
bodeplot(Ps, w, legend=false)
```
We can also convert the vector of system models to a single model with [`RobustAndOptimalControl.ss2particles`](@ref), which will convert the coefficients of the state space models to [`MonteCarloMeasurements.Particles`](https://baggepinnen.github.io/MonteCarloMeasurements.jl/latest/) objects.
```@example BATCHLIN
using RobustAndOptimalControl
P = RobustAndOptimalControl.ss2particles(Ps) # convert to a single StateSpace system with `Particles` as coefficients.
```

notice how some coefficients are plotted like uncertain numbers `-19.2 ± 7.6`. We can plot such models as well:
```@example BATCHLIN
bodeplot(P, w, legend=:bottomright) # Should look similar to the one above
```

## Controller tuning
Let's also do some controller tuning for the linearized models above. The function `batch_tune` is not really required here, but it shows how we might go about building more sophisticated tools for batch tuning. In this example, we will tune a PID controller using the function `loopshapingPID`.
```@example BATCHLIN
function batch_tune(f, Ps)
    f.(Ps)
end

Cs = batch_tune(Ps) do P
    C, kp, ki, kd, fig, CF = loopshapingPID(P, 7; Mt=1.2, Tf = 1/100)
    ss(CF)
end

P = RobustAndOptimalControl.ss2particles(Ps)
C = RobustAndOptimalControl.ss2particles(Cs)

nyquistplot(P * C,
            w,
            ylims = (-10, 2),
            xlims = (-8, 5),
            points = true,
            Ms_circles = [1.5, 2],
            Mt_circles = [1.2])
```
Above, we plotted the Nyquist curve of the loop-transfer function for all system realizations. RobustAndOptimalControl.jl has some facilities for fitting circles around the Nyquist curve for uncertain systems, which we could use here:
```@example BATCHLIN
centers, radii = fit_complex_perturbations(P * C, w; relative = false, nominal = :center)
nyquistcircles!(w, centers, radii, ylims = (-5, 1), xlims = (-3, 4))
```
some methods for robust control operate on such circles. Notice how the circles are conservative in many cases, this is typically due to the gain varying between the models for the same phase.

If you plot the Nyquist curve using the `plotly()` backend rather than the default `gr()` backend used here, you can hover the mouse over the curves and see which frequency they correspond to etc. 

## Gain scheduling
Above, we tuned one controller for each operating point, wouldn't it be nice if we had some features to simulate a [gain-scheduled controller](https://en.wikipedia.org/wiki/Gain_scheduling) that interpolates between the different controllers depending on the operating pont? [`GainScheduledStateSpace`](@ref) is such a thing, we show how to use it below. For fun, we simulate some reference step responses for each individual controller in the array `Cs` and end with simulating the gain-scheduled controller.

```@example BATCHLIN
using OrdinaryDiffEq
using DataInterpolations # Required to interpolate between the controllers
@named fb  = Blocks.Add(k2=-1)
@named ref = Blocks.Square(frequency=1/6, amplitude=0.5, offset=0.5, start_time=1)
@named F   = Blocks.SecondOrder(w=15, d=1) # A reference pre-filter
connect    = ModelingToolkit.connect

closed_loop_eqs = [
    connect(ref.output, F.input)
    connect(F.output, fb.input1)
    connect(duffing.y, fb.input2)
]
plot(layout=2)

# Simulate each individual controller
for C in Cs
    @named Ci = ODESystem(C)
    eqs = [
        closed_loop_eqs
        connect(fb.output, Ci.input)
        connect(Ci.output, duffing.u)
    ]
    @named closed_loop = ODESystem(eqs, t, systems=[duffing, Ci, fb, ref, F])
    prob = ODEProblem(structural_simplify(closed_loop), [], (0.0, 8.0))
    sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
    plot!(sol, idxs=[duffing.y.u, duffing.u.u], layout=2, lab="")
end

# Simulate gain-scheduled controller
@named Cgs = GainScheduledStateSpace(Cs, xs, interpolator=LinearInterpolation)
eqs = [
    closed_loop_eqs
    connect(fb.output, Cgs.input)
    connect(Cgs.output, duffing.u)
    connect(duffing.y, Cgs.scheduling_input) # Don't forget to connect the scheduling variable!
]
@named closed_loop = ODESystem(eqs, t, systems=[duffing, Cgs, fb, ref, F])
prob = ODEProblem(structural_simplify(closed_loop), [], (0.0, 8.0))
sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
plot!(sol, idxs=[duffing.y.u, duffing.u.u], l=(2, :red), lab="Gain scheduled")
plot!(sol, idxs=F.output.u, l=(1, :black, :dash, 0.5), lab="Ref")
```

If everything worked as expected, the gain-scheduled controller should perform better than each of the included controllers individually. 


## Summary
We have seen how to
- Perform linearization of a nonlinear ModelingToolkit model in multiple different operating points
- Handle arrays of models or models with `Particles` as coefficients
- Simulate a gain-scheduled controller that interpolates between linear controllers


Batch linearization in multiple different operating points is an intuitive way to perform analysis of a nonlinear control system. Gain-scheduling is an equally intuitive way of realizing a nonlinear controller. Care should be taken to ensure that the scheduling variable does not change too fast such that the linear assumption at each instance of time is not violated.


```@example BATCHLIN
using Test
@test sol(6.99, idxs=closed_loop.duffing.y.u) ≈ 0.0 atol=0.01
```