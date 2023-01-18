# Batch Linearization and gain scheduling
This example will demonstrate how to linearize a nonlinear ModelingToolkit model in multiple different operating points, and some tools to work with groups of linear models representing the same system in different operating points.


## System model
The model will be a simple Duffing oscillator:
```@example BATCHLIN
using ControlSystemsMTK, ModelingToolkit, MonteCarloMeasurements
using ModelingToolkit: getdefault
unsafe_comparisons(true)

# Create a model
@parameters t k=10 k3=2 c=1
@variables x(t)=0 [bounds = (-2, 2)]
@variables v(t)=0
@variables u(t)=0
@variables y(t)=0

D = Differential(t)

eqs = [D(x) ~ v
       D(v) ~ -k * x - k3 * x^3 - c * v + 10u
       y ~ x]


@named duffing = ODESystem(eqs, t)
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
Ps, ssys = batch_ss(duffing, [u], [y], ops)
nothing # hide
```

Plotting functions like [`bodeplot`](@ref) accept vectors of systems, so this works
```@example BATCHLIN
using ControlSystemsBase, Plots
w = exp10.(LinRange(-2, 3, 200))
bodeplot(Ps, w)
```
We can also convert the vector of system models to a single model with [`RobustAndOptimalControl.ss2particles`](@ref), which will convert the coefficients of the state space models to [`MonteCarloMeasurements.Particles`](https://baggepinnen.github.io/MonteCarloMeasurements.jl/latest/) objects.
```@example BATCHLIN
using RobustAndOptimalControl
P = RobustAndOptimalControl.ss2particles(Ps) # convert to a single StateSpace system with `Particles` as coefficients.
```

notice how some coefficients are plotted like uncertain numbers `-19.2 ± 7.6`. We can plot such models as well:
```@example BATCHLIN
bodeplot(P, w) # Should look similar to the one above
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
Above, we tuned one controller for each operating point, wouldn't it be nice if we had some features to simulate a gain-scheduled controller that interpolates between the different controllers depending on the operating pont? Currently, we don't have that, so maybe we have to mention this to Santa next time we see him. 