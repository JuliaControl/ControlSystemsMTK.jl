# Batch Linearization and gain scheduling
This example will demonstrate how to linearize a nonlinear ModelingToolkit model in multiple different operating points, and some tools to work with groups of linear models representing the same system in different operating points. We'll end with designing and simulating a [gain-scheduled controller](https://en.wikipedia.org/wiki/Gain_scheduling), i.e., a nonlinear controller created as an interpolation between linear controllers.

!!! note "What is an operating point?"
    An operating point is typically understood as a tuple of the form ``(x, u)``, where ``x`` is the state vector and ``u`` is the input vector. However, we may choose to include *parameters* ``p`` in the operating point as well. This may be useful when some parameters are uncertain or time varying, and we want to perform analysis over multiple possible parameter values.

## System model
The model will be a simple Duffing oscillator:
```@example BATCHLIN
using ControlSystemsMTK, ModelingToolkit, MonteCarloMeasurements, ModelingToolkitStandardLibrary.Blocks
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


@named duffing = System(eqs, t, systems=[y, u])
```

## Batch linearization
To perform batch linearization, we create a vector of operating points, and then linearize the model around each of these points. The function [`batch_ss`](@ref) does this for us, and returns a vector of `StateSpace` models, one for each operating point. An operating point is a `Dict` that maps variables in the MTK model to numerical values. In the example below, we simply sample the variables uniformly within their bounds specified when we created the variables (normally, we might want to linearize on stationary points)
```@example BATCHLIN
N = 5 # Number of samples
xs = range(getbounds(x)[1], getbounds(x)[2], length=N)
ops = Dict.(u.u => 0, x .=> xs)
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

notice how some coefficients are plotted like uncertain numbers `-13.8 ± 4.3`. We can plot such models as well:
```@example BATCHLIN
bodeplot(P, w, legend=:bottomright, adaptive=false) # Should look similar to the one above
```

## Controller tuning
Let's also do some controller tuning for the linearized models above. The function `batch_tune` is not really required here, but it shows how we might go about building more sophisticated tools for batch tuning. In this example, we will tune a PID controller using the function [`loopshapingPID`](@ref). Note, this procedure is not limited to tuning a gain-scheduled PID controller, it should work for gain-scheduling of any LTI controller. 

!!! "note" Interpolating between controllers
    There are multiple ways in which one could interpolate between different controllers. The two most common approaches are to interpolate their outputs, and interpolating their coefficients. When interpolating the coefficients, it is important to ensure that all controllers have the same structure for the interpolation to be meaningful. One may for example interpolate between PID coefficients, or between the coefficients of a state-space model. When interpolating state-space matrices, the systems must all share the same basis, i.e., the state variables must all have the same meaning among the interpolated systems. When converting a transfer function to state-space form, a numerical balancing is performed, this alters the meaning of the state variables and may introduce artificial dynamics to the interpolated system. We thus pass `balance=false` to the function `ss` to avoid this, or pick a form explicitly, e.g., `modal_form`.

```@example BATCHLIN
function batch_tune(f, Ps)
    f.(Ps)
end

Cs = batch_tune(Ps) do P
    C, kp, ki, kd, fig, CF = loopshapingPID(P, 7; Mt=1.2, Tf = 1/100)
    ss(CF, balance=false)
    modal_form(ss(CF, balance=true))[1]
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
using OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqRosenbrock
using DataInterpolations # Required to interpolate between the controllers
@named fb  = Blocks.Add(k2=-1)
@named ref = Blocks.Square(frequency=1/6, amplitude=0.5, offset=0.5, start_time=1)
@named F   = Blocks.SecondOrder(w=15, d=1) # A reference pre-filter
connect    = ModelingToolkit.connect

closed_loop_eqs = [
    connect(ref.output, :r0, F.input)
    connect(F.output,  :r, fb.input1) # Add an analysis point :r
    connect(duffing.y, :y, fb.input2) # Add an analysis point :y
]
plot(layout=2)

# Simulate each individual controller
for C in Cs
    @named Ci = System(C)
    eqs = [
        closed_loop_eqs
        connect(fb.output, Ci.input)
        connect(Ci.output, duffing.u)
    ]
    @named closed_loop = System(eqs, t, systems=[duffing, Ci, fb, ref, F])
    prob = ODEProblem(structural_simplify(closed_loop), [F.x => 0, F.xd => 0], (0.0, 8.0))
    sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
    plot!(sol, idxs=[duffing.y.u, duffing.u.u], layout=2, lab="")
end

# Simulate gain-scheduled controller
@named Cgs = GainScheduledStateSpace(Cs, xs, interpolator=LinearInterpolation)
eqs = [
    closed_loop_eqs
    connect(fb.output, Cgs.input)
    connect(Cgs.output, :u, duffing.u)
    connect(duffing.y, :v, Cgs.scheduling_input) # Don't forget to connect the scheduling variable!
]
@named closed_loop = System(eqs, t, systems=[duffing, Cgs, fb, ref, F])
prob = ODEProblem(structural_simplify(closed_loop), [F.xd => 0], (0.0, 8.0))
sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8, initializealg=SciMLBase.NoInit(), dtmax=0.01)
plot!(sol, idxs=[duffing.y.u, duffing.u.u], l=(2, :red), lab="Gain scheduled")
plot!(sol, idxs=F.output.u, l=(1, :black, :dash, 0.5), lab="Ref")
```

If everything worked as expected, the gain-scheduled controller should perform reasonably well across the entire scheduling range. 


## C-Code generation
We can generate C-code to interpolate our controller using the function [`SymbolicControlSystems.print_c_array`](@ref) from [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl). If the controller is a standard [`ControlSystemsBase.StateSpace`](@ref) object, a function that filters the input through the controller can be generated by calling [`SymbolicControlSystems.ccode`](@ref). But if the controller is a vector of controllers representing a gain-scheduled controller, a function that creates the interpolated dynamics is written. In the code below, we shorten the vector of controllers to make the generated C-code easier to read by passing `Cs[1:3:end]` and `xs[1:3:end]`
```@example BATCHLIN
using SymbolicControlSystems, ControlSystemsBase
Cs_disc = c2d.(Cs, 0.05, :tustin) # Discretize the controller before generating code
code = SymbolicControlSystems.print_c_array(stdout, Cs_disc[1:3:end], xs[1:3:end], "Cs")
```
The generated code starts by defining the interpolation vector `xs`, this variable is called `Cs_interp_vect` in the generated code. The code then defines all the ``A`` matrices as a 3-dimensional array, followed by a function that performs the interpolation `interpolate_Cs_A`. This function takes the output array as the first argument, a pointer to the 3D array with interpolation matrices, the interpolation vector as well as the interpolation variable `t`, in this document called ``v``. The same code is then repeated for the matrices ``B,C,D`` as well if they require interpolation (if they are all the same, no interpolation code is written). 

## Linearize around a trajectory
We can linearize around a trajectory obtained from `solve` using the function [`trajectory_ss`](@ref). We provide it with a vector of time points along the trajectory at which to linearize, and in this case we specify the inputs and outputs to linearize between as analysis points `r` and `y`.
```@example BATCHLIN
timepoints = 0:0.01:8
Ps2, ssys, ops2, resolved_ops = trajectory_ss(closed_loop, closed_loop.r, closed_loop.y, sol; t=timepoints, verbose=true);
bodeplot(Ps2, w, legend=false)
```
Not how the closed-loop system changes very little along the trajectory, this is a good indication that the gain-scheduled controller is able to make the system appear linear.

Internally, [`trajectory_ss`](@ref) works very much the same as [`batch_ss`](@ref), but constructs operating points automatically along the trajectory. This requires that the solution contains the states of the simplified system, accessible through the `idxs` argument like `sol(t, idxs=x)`. By linearizing the same system as we simulated, we ensure that this condition holds, doing so requires that we specify the inputs and outputs as analysis points rather than as variables.


We can replicate the figure above by linearizing the plant and the controller individually, by providing the `loop_openings` argument. When linearizing the plant, we disconnect the controller input by passing `loop_openings=[closed_loop.u]`, and when linearizing the controller, we have various options for disconnecting the the plant:
- Break the connection from plant output to controller input by passing `loop_openings=[closed_loop.y]`
- Break the connection between the controller and the plant input by passing `loop_openings=[closed_loop.u]`
- Break the connection `y` as well as the scheduling variable `v` (which is another form of feedback) by passing `loop_openings=[closed_loop.y, closed_loop.v]`

We will explore these options below, starting with the first option, breaking the connection `y`:
```@example BATCHLIN
kwargs = (; adaptive=false, legend=false)
plants, _ = trajectory_ss(closed_loop, closed_loop.u, closed_loop.y, sol; t=timepoints, verbose=true, loop_openings=[closed_loop.u]);
controllersy, ssy, ops3, resolved_ops3 = trajectory_ss(closed_loop, closed_loop.r, closed_loop.u, sol; t=timepoints, verbose=true, loop_openings=[closed_loop.y]);

closed_loopsy = feedback.(plants .* controllersy)
bodeplot(closed_loopsy, w; title="Loop open at y", kwargs...)
```
When we open the loop at `u`, we get a slightly different result:
```@example BATCHLIN
controllersu, ssu = trajectory_ss(closed_loop, closed_loop.r, closed_loop.u, sol; t=timepoints, verbose=true, loop_openings=[closed_loop.u]);

closed_loopsu = feedback.(plants .* controllersu)
bodeplot(closed_loopsu, w; title="Loop open at u", kwargs...)
```
In this case, all static gains are 1. A similar result is obtained by explicitly breaking the scheduling feedback `v` in addition to an opening of `y`:
```@example BATCHLIN
controllersv, ssv = trajectory_ss(closed_loop, closed_loop.r, closed_loop.u, sol; t=timepoints, verbose=true, loop_openings=[closed_loop.y, closed_loop.v]);

closed_loopsv = feedback.(plants .* controllersv)
bodeplot(closed_loopsv, w; title="Loop open at v and y", kwargs...)
```

We have thus far treated the controller as a SISO system, but we could also view it as a system with two inputs, the measurement feedback and the scheduling feedback. For completeness, we show below how to derive the corresponding MIMO systems

```@example BATCHLIN
plants_mimo, _ = trajectory_ss(closed_loop, closed_loop.u, [closed_loop.y, closed_loop.v], sol; t=timepoints, verbose=true, loop_openings=[closed_loop.u]);
controllers_mimo, ssm = trajectory_ss(closed_loop, [closed_loop.r, closed_loop.v], closed_loop.u, sol; t=timepoints, verbose=true, loop_openings=[closed_loop.u]);

closed_loops_mimo = feedback.(controllers_mimo .* plants_mimo) # Look at complementary sensitivity function in the input, since this is a SISO system
bodeplot(closed_loops_mimo, w; title="Loop open at MIMO", kwargs...)
```



Why are the results different depending on where we open the loop? We can understand the difference by comparing the Bode plots of the controllers. 

```@example BATCHLIN
plot(
    bodeplot(Cs,           w, legend=false, plotphase=false, title="Designed controllers"),
    bodeplot(controllersy, w, legend=false, plotphase=false, title="Loop open at y"),
    bodeplot(controllersu, w, legend=false, plotphase=false, title="Loop open at u"),
    bodeplot(controllersv, w, legend=false, plotphase=false, title="Loop open at v and y"),
)
```
if we open at both `y` and `v` or we open at `u`, we get controllers for the different values of the scheduling variable, and the corresponding measurement feedback (which is the same as the scheduling variable in this case).
```@example BATCHLIN
using Test
@test all(sminreal.(controllersv) .== sminreal.(controllersu))
```

However, if we only open at `y` we get controller linearizations that _still contain the closed loop through the scheduling connection_ `v`. We can verify this by looking at what variables are present in the input-output map
```@example BATCHLIN
sminreal(controllersy[end]).x
```
notice how the state of the plant is included in the controller, this is an indication that we didn't fully isolate the controller during the linearizaiton. If we do the same thing for the controller with the loop opened at `u`, we see that the state of the plant is not included in the controller:
```@example BATCHLIN
sminreal(controllersu[end]).x
```
The call to `sminreal` is important here, it removes the states that are not needed to represent the input-output map of the system. The state of the full model, including the plant state, is present in the linearization before this call. 



The easiest way to ensure that the controller is properly disconnected from the plant while taking into account the different scheduling along the trajectory is thus to break at `u`.

## Summary
We have seen how to
- Perform linearization of a nonlinear ModelingToolkit model in multiple different operating points
- Handle arrays of models or models with `Particles` as coefficients
- Simulate a gain-scheduled controller that interpolates between linear controllers
- Write C-code to perform the interpolation of the controllers


Batch linearization in multiple different operating points is an intuitive way to perform analysis of a nonlinear control system. Gain-scheduling is an equally intuitive way of realizing a nonlinear controller. Care should be taken to ensure that the scheduling variable does not change too fast such that the linear assumption at each instance of time is not violated.


```@example BATCHLIN
@test sol(6.99, idxs=closed_loop.duffing.y.u) ≈ 0.0 atol=0.01
```