# Hybrid and multi-rate models

A *hybrid* model contains clocked (discrete-time) partitions, such as digital controllers, that are connected to continuous-time equations through the operators `Sample` and `Hold`. This page describes how such models are linearized and analyzed:

- With a single clock, the model sampled at the clock is time invariant, and [`named_ss`](@ref) and the `get_named_*` functions return a discrete-time system.
- With several clocks, the model is *periodically time varying*. It has no transfer function, and the functions on this page provide the representations and the frequency-domain quantities that replace it.

The linearization of each clock partition is performed by `ModelingToolkit.linearize_hybrid`.

## A model with one clock

The model below consists of a continuous-time plant and a clocked PI controller. The controller output is held and connected to the plant through the analysis point `plant_input`.

```@example MULTIRATE
using ControlSystemsMTK, ControlSystemsBase, ModelingToolkit, RobustAndOptimalControl, LinearAlgebra, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D, connect
default(size = (700, 400))

Tfast = 0.125          # Period of the fast clock
ratio = 4
Tslow = ratio * Tfast  # Period of the slow clock
kslow = ShiftIndex(Clock(Tslow))
kfast = ShiftIndex(Clock(Tfast))

@variables x1(t) x2(t) y(t) d(t) r(t) e(t) xi(t) ud(t) yf(t)
@variables uh(t) [output = true] u(t) [input = true]
@parameters kp = 0.4 ki = 0.3

plant = [
    D(x1) ~ x2
    D(x2) ~ -x1 - 0.6x2 + u + d
    y ~ x1
]
# A PI controller on the slow clock, `measurement` is a signal on that clock
controller(measurement) = [
    e ~ r - measurement
    xi(kslow) ~ xi(kslow - 1) + ki * Tslow * e(kslow)
    ud ~ kp * e + xi
]
hold = [uh ~ Hold(ud); connect(uh, :plant_input, u)]

@named single = System([plant; controller(Sample(Tslow)(y)); hold; r ~ 0], t)
op = Dict(x1 => 0.0, x2 => 0.0, d => 0.0, xi => 0.0, r => 0.0, yf => 0.0)
G = named_ss(single, [d], [y]; op, hybrid = true)
```

The result is the system from the disturbance `d` to the output `y` with the loop closed: the continuous-time partition is discretized with zero-order hold at the period of the clock and connected to the linearization of the controller. The disturbance is consequently assumed to be constant between the ticks of the clock, and the output is obtained at the ticks only. The keyword argument `hybrid = true` selects this linearization. It is the default when a discrete compilation pass such as `SynchToolkit.compile_lustre` is supplied in `additional_passes`, which is how hybrid models are compiled for simulation.

The sensitivity functions are obtained in the same way, here for a model in which the disturbance is absent:

```@example MULTIRATE
@named single_closed = System([plant; controller(Sample(Tslow)(y)); hold; r ~ 0; d ~ 0], t)
S = get_named_sensitivity(single_closed, single_closed.plant_input; op, hybrid = true)
bodeplot(S, title = "Sensitivity at the plant input")
```

## A model with two clocks

The measurement is now filtered by a filter on a fast clock, and the controller on the slow clock reads the output of the filter.

```@example MULTIRATE
fastfilter = [yf(kfast) ~ 0.5 * yf(kfast - 1) + 0.5 * Sample(Tfast)(y)]
slow(signal) = Sample(Tslow)(Hold(signal))

@named multi = System([plant; fastfilter; controller(slow(yf)); hold], t)
nothing # hide
```

`named_ss` rejects this model, since no time-invariant system represents it:

```@example MULTIRATE
try
    named_ss(multi, [d, r], [y]; op, hybrid = true)
catch err
    println(err.msg)
end
```

### Why the model is periodic

The slow partition is updated at every fourth tick of the fast clock. The response of the model to an input therefore depends on the time at which the input is applied in relation to the ticks of the slow clock, but a shift of the input by the common period `Tc` of all clocks shifts the response by `Tc`. A model with this property is periodically time varying with period `Tc`.

The consequence in the frequency domain is that a sinusoidal input with angular frequency ``\omega`` gives a steady-state output that contains the frequencies

```math
\omega + k\,\dfrac{2\pi}{T_c}, \qquad k = 0, \pm 1, \pm 2, \ldots
```

The component with ``k = 0`` is the *fundamental* component, the others are *harmonics* created by sampling and holding (aliasing). A time-invariant system produces the fundamental component only.

### The periodic model and the lifted system

[`multirate_linearization`](@ref) creates a linear periodic model on the *event grid*, the sorted set of instants within one period at which a clock ticks or a signal of the continuous-time partition is updated. The period is the least common multiple of the clock periods, which are required to be commensurate.

```@example MULTIRATE
ml = multirate_linearization(multi, [d, r], [y, ud]; op)
```

The model is a sequence of maps, one per event ``k``,

```math
\begin{aligned}
x_{k+1} &= A_k x_k + B_k u \\
y_k &= C_k x_k + D_k u
\end{aligned}
```

where the state ``x`` contains the state of the continuous-time partition, the states of the clocked partitions and the held values of the signals that cross clock boundaries. At an instant where clocks tick:

1. A sampled signal of the continuous-time partition is read before any hold is updated at that instant.
2. A clocked partition that reads a signal of another clocked partition that ticks at the same instant reads the new value.
3. All holds are updated, and the continuous-time partition is integrated exactly to the next event with constant inputs.

Each signal has a set of instants within the period: the ticks of its clock for a signal of a clocked partition, and a uniform grid for a signal of the continuous-time partition. An input is held between its instants, an output is observed at them.

*Lifting* collects the values of every signal over one period into a vector. The periodic model then becomes the time-invariant system

```math
\begin{aligned}
x(m+1) &= F x(m) + G \bar u(m) \\
\bar y(m) &= H x(m) + E \bar u(m)
\end{aligned}
```

with sample interval `Tc`, where ``\bar u(m)`` and ``\bar y(m)`` contain all values of the signals in period ``m``, and

```math
F = A_K \cdots A_2 A_1
```

is the *monodromy matrix*.

```@example MULTIRATE
Gl = lifted_ss(ml)
Gl.u, Gl.y
```

The disturbance `d` and the output `y` belong to the continuous-time partition and have one channel per tick of the fast clock, the reference `r` and the controller output `ud` have one channel per period. The lifted system is an ordinary discrete-time system and can be simulated, for example for a step in the reference:

```@example MULTIRATE
nperiods = 40
ubar = zeros(Gl.nu, nperiods)
ubar[end, :] .= 1 # The reference is the last channel
res = lsim(Gl.sys, ubar)
tfast = Tfast .* (0:ratio*nperiods-1)
plot(tfast, vec(res.y[1:ratio, :]), label = "y", xlabel = "Time")
plot!(Tslow .* (0:nperiods-1), res.y[end, :], label = "ud", linetype = :steppost)
```

### Stability

The eigenvalues of the monodromy matrix are the *characteristic multipliers* of the periodic model. The model is exponentially stable if and only if all multipliers have a magnitude smaller than one. This is the only stability test on this page that is exact without qualification.

```@example MULTIRATE
λ = multipliers(ml)
maximum(abs, λ)
```

The hold variables contribute multipliers that are zero.

### Frequency-domain functions

For an input port ``u`` with instants ``t_{u,1}, \ldots, t_{u,n_u}`` and an output port ``y`` with instants ``t_{y,1}, \ldots, t_{y,n_y}``, let ``G_L`` be the corresponding block of the transfer matrix of the lifted system and

```math
a_s(\omega) = \begin{bmatrix} e^{i\omega t_{s,1}} & \cdots & e^{i\omega t_{s,n_s}} \end{bmatrix}^{\mathsf T}
```

the values of the sinusoid ``e^{i\omega t}`` at the instants of a port ``s``. The following functions of ``\omega`` are available.

| Function | Definition | Interpretation |
|:---|:---|:---|
| [`fundamental_response`](@ref) | ``a_y(\omega)^* G_L(e^{i\omega T_c})\, a_u(\omega) / n_y`` | Complex amplitude of the fundamental component of the output |
| [`harmonic_response`](@ref) | ``a_y(\omega + 2\pi k/T_c)^* G_L(e^{i\omega T_c})\, a_u(\omega) / n_y`` | Complex amplitude of harmonic ``k`` |
| [`performance_frequency_gain`](@ref) | ``\lVert G_L(e^{i\omega T_c})\, a_u(\omega)\rVert_2 / \sqrt{n_y}`` | Root-mean-square gain for one sinusoid, all harmonics included |
| [`robust_frequency_gain`](@ref) | ``\bar\sigma\big(G_L(e^{i\omega T_c})\big)`` | Largest gain over inputs that contain several of the frequencies ``\omega + 2\pi k/T_c`` |

They satisfy

```math
|\text{fundamental response}| \le \text{performance frequency gain} \le \text{robust frequency gain}
```

for ports with the same number of instants, and the squared performance frequency gain is the sum of the squared magnitudes of all harmonics. For a time-invariant system the three quantities are equal to the magnitude of the frequency response.

```@example MULTIRATE
w = exp10.(range(-1, log10(π / Tfast), length = 400))
plot(w, abs.(fundamental_response(ml, "d(t)", "y(t)", w)), label = "Fundamental response",
    xscale = :log10, yscale = :log10, xlabel = "Frequency [rad/s]", title = "From d to y", legend = :bottomleft)
plot!(w, performance_frequency_gain(ml, "d(t)", "y(t)", w), label = "Performance frequency gain")
plot!(w, robust_frequency_gain(ml, "d(t)", "y(t)", w), label = "Robust frequency gain")
vline!([π / Tslow], label = "Nyquist frequency of the slow clock", linestyle = :dash, color = :black)
```

The fundamental response and the performance frequency gain agree at frequencies where the harmonics are small. The ratio between them is a measure of the amount of aliasing, and a ratio close to one indicates that the model behaves as a time-invariant system for a sinusoid of that frequency. The harmonics themselves show to which frequencies the input is converted:

```@example MULTIRATE
plot(xscale = :log10, yscale = :log10, xlabel = "Frequency of the input [rad/s]", title = "Harmonics of y for a sinusoidal d", legend = :bottomleft, ylims = (1e-6, 10))
for k in 0:ratio-1
    plot!(w, abs.(harmonic_response(ml, "d(t)", "y(t)", w, k)), label = "k = $k")
end
current()
```

The following properties determine how the functions may be used.

- The fundamental response is not a transfer function: it is not the ratio between the transforms of the output and the input. Neither stability nor stability margins can be concluded from it. For the sensitivity functions of a model with one clock, it is the *fundamental sensitivity function* of Goodwin and Salgado (1994), for which a large peak is a proof of poor robustness, whereas a small peak is not a proof of good robustness.
- The supremum over frequency of the robust frequency gain is the induced norm between the sequences of sampled signals. It is periodic in ``\omega`` with period ``2\pi/T_c``.
- The functions describe the signals at their instants. The behavior of a continuous-time signal between the ticks of the clocks is made visible by a finer grid, see the next section.

### Behavior between the ticks

Signals of the continuous-time partition are represented on a uniform grid, the interval of which is the smallest clock period divided by the keyword argument `oversample`. A larger value refines the grid without changing the ticks of the clocks, such that the multipliers are unchanged while the frequency functions extend to higher frequencies and include the response between the ticks. `continuous_rate = :period` instead selects the common period as interval, which gives one channel per signal.

```@example MULTIRATE
ml4 = multirate_linearization(multi, [d, r], [y, ud]; op, oversample = 4)
w4 = exp10.(range(-1, log10(4π / Tfast), length = 400))
plot(w, performance_frequency_gain(ml, "d(t)", "y(t)", w), label = "oversample = 1",
    xscale = :log10, yscale = :log10, xlabel = "Frequency [rad/s]", title = "Performance frequency gain from d to y", legend = :bottomleft)
plot!(w4, performance_frequency_gain(ml4, "d(t)", "y(t)", w4), label = "oversample = 4")
```

```@example MULTIRATE
maximum(abs, multipliers(ml4)) ≈ maximum(abs, multipliers(ml))
```

As the grid is refined, the frequency functions approach those of the continuous-time signals. No bound on the difference is available for a finite `oversample`.

### Sensitivity functions

`get_sensitivity`, `get_comp_sensitivity` and `get_looptransfer` return a `HybridLinearization` with the keyword argument `hybrid = true`, from which the periodic model is created. The sensitivity function at the plant input is

```@example MULTIRATE
@named multi_closed = System([plant; fastfilter; controller(slow(yf)); hold; r ~ 0; d ~ 0], t)
hl = get_sensitivity(multi_closed, multi_closed.plant_input; op, hybrid = true)
mlS = multirate_linearization(hl)
plot(w, abs.(fundamental_response(mlS, 1, 1, w)), label = "Fundamental sensitivity",
    xscale = :log10, yscale = :log10, xlabel = "Frequency [rad/s]", legend = :bottomright)
plot!(w, performance_frequency_gain(mlS, 1, 1, w), label = "Performance frequency gain")
plot!(w, robust_frequency_gain(mlS, 1, 1, w), label = "Robust frequency gain")
```

### Stability margins

In a time-invariant loop, a gain or phase margin states how much one perturbation of the loop gain is tolerated. In a periodic loop, the signal at the break point has ``n`` instants per period, and perturbations differ in how they act on these ``n`` values. [`multirate_margins`](@ref) returns one margin per class of perturbations:

| Field | Perturbation | Property |
|:---|:---|:---|
| `gain_interval` | A constant real gain | Exact interval of gains for which the loop is stable |
| `disk_margin_lti_lower_bound` | A time-invariant dynamic system at the rate of the signal | A lower bound of the disk margin on a frequency grid |
| `disk_margin_dynamic` | An arbitrary dynamic system that acts on the lifted signal, including periodically time-varying systems | Exact for this class, conservative for time-invariant perturbations |

A time-invariant perturbation ``\Delta`` at the rate of the signal acts on the frequencies ``\omega + 2\pi k/T_c``, ``k = 0, \ldots, n-1``, with the values ``\Delta(e^{i(\omega + 2\pi k/T_c)T_c/n})``, which are independent of each other. In the basis of these frequencies the perturbation is diagonal, and the margin follows from the structured singular value with ``n`` independent complex scalars, for which an upper bound is computed. The disk margins are `RobustAndOptimalControl.Diskmargin` objects, see [the documentation of disk margins](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#Diskmargin-example).

The example below places an analysis point at the output of the filter, a signal on the fast clock with four instants per period.

```@example MULTIRATE
@variables fo(t) [output = true] fi(t) [input = true] yf2(t)
@named multi_break = System([plant; fastfilter; fo ~ yf; connect(fo, :filter_output, fi);
    yf2(kfast) ~ fi(kfast); controller(slow(yf2)); hold; r ~ 0; d ~ 0], t)
opb = merge(op, Dict(yf2 => 0.0))
m = multirate_margins(multi_break, multi_break.filter_output; op = opb)
m.gain_interval
```

```@example MULTIRATE
m.disk_margin_lti_lower_bound
```

```@example MULTIRATE
m.disk_margin_dynamic
```

The larger the class of perturbations, the smaller the margin. If the break point has one instant per period, the loop seen from the break point is time invariant, and the two disk margins coincide with `RobustAndOptimalControl.diskmargin` of the loop-transfer function.

For an analysis point in the continuous-time partition, the perturbations act on the values of the signal on its grid. Only the gain interval then refers to a perturbation of the continuous-time signal itself; margins for dynamic perturbations of a continuous-time signal, such as a delay margin, are not provided.

### A time-invariant approximation

[`single_rate_ss`](@ref) converts all partitions to one common rate and connects them, which gives an ordinary state-space system at the cost of an approximation: aliasing between the rates is not represented. This is the approach of block-wise rate conversion that is common in simulation tools. With `Ts = 0` the result is a continuous-time system in which the clocked partitions are converted with the Tustin approximation.

```@example MULTIRATE
Gc = single_rate_ss(multi, [d, r], [y, ud]; Ts = 0, op)
approx = abs.(vec(freqresp(Gc.sys[1, 1], w)))
plot(w, abs.(fundamental_response(ml4, "d(t)", "y(t)", w)), label = "Fundamental response",
    xscale = :log10, yscale = :log10, xlabel = "Frequency [rad/s]", title = "From d to y", legend = :bottomleft)
plot!(w, approx, label = "single_rate_ss, continuous time", linestyle = :dash)
vline!([π / Tslow], label = "Nyquist frequency of the slow clock", linestyle = :dash, color = :black)
```

The approximation agrees with the fundamental response at low frequencies. It deviates around the resonance of the closed loop, because the delay that the hold introduces is not represented in the continuous-time result, and it has no validity above the Nyquist frequency of the slowest clock. No stability or robustness property of the model follows from it; in particular, the poles of the approximation are not the multipliers of the model. It is the only representation on this page that is available when the clock periods are not commensurate.

## Choice of representation

| Task | Function | Validity |
|:---|:---|:---|
| Stability | [`multipliers`](@ref) | Exact |
| Simulation, norms, synthesis with tools for time-invariant systems | [`lifted_ss`](@ref) | Exact at the instants of the signals |
| Loop shaping, location of a peak in frequency | [`fundamental_response`](@ref), [`harmonic_response`](@ref) | Exact amplitude of each component; no certificate of stability or robustness |
| Detection of aliasing | [`performance_frequency_gain`](@ref) in relation to the fundamental response | Exact at the instants of the signals |
| Largest gain | [`robust_frequency_gain`](@ref) | Induced norm between the sampled sequences |
| Margins | [`multirate_margins`](@ref) | Per class of perturbations, see above |
| A rational model at one rate | [`single_rate_ss`](@ref) | Approximation |

## Limitations

- The linearization is required to be performed in a stationary operating point of the whole model: an equilibrium of the continuous-time partition at which the states of the clocked partitions are constant. This is not verified. A model that settles on a periodic orbit rather than an equilibrium is not supported.
- The clocks are required to be periodic with commensurate periods. Events, state machines and clocks without a fixed period are not accounted for.
- Ticks of different clocks that coincide are treated as simultaneous. In a simulation, tick times computed in floating-point arithmetic coincide only if the periods are exactly representable (for example 0.125 and 0.5, as opposed to 0.1 and 0.4), and partitions are otherwise evaluated in the order of their tick times.
- The connection between two clocked partitions is expressed as `Sample(Hold(x))` in ModelingToolkit and as `Latest` in SynchToolkit. Both are treated according to item 2 of the semantics above, which was verified against simulations with SynchToolkit.
- With analysis points, all inputs, and all outputs, are required to belong to one clock partition.

## References

- G. C. Goodwin and M. Salgado, "Frequency domain sensitivity functions for continuous time systems under sampled data control", *Automatica* 30(8), 1994.
- M. M. Seron, J. H. Braslavsky and G. C. Goodwin, *Fundamental Limitations in Filtering and Control*, Springer, 1997, chapter 6.
- R. A. Meyer and C. S. Burrus, "A unified analysis of multirate and periodically time-varying digital filters", *IEEE Transactions on Circuits and Systems* 22(3), 1975.
- M. Araki and K. Yamamoto, "Multivariable multirate sampled-data systems: state-space description, transfer characteristics, and Nyquist criterion", *IEEE Transactions on Automatic Control* 31(2), 1986.
- O. Lindgärde and B. Lennartson, "Performance and robust frequency response for multirate sampled-data systems", *American Control Conference*, 1997.
- M. van Haren, L. Blanken and T. Oomen, "Performance analysis of multirate systems: a direct frequency-domain identification approach", *Mechanical Systems and Signal Processing* 235, 2025.
- G. Dullerud and K. Glover, "Robust stabilization of sampled-data systems to structured LTI perturbations", *IEEE Transactions on Automatic Control* 38(10), 1993.
- P. M. Thompson, "Gain and phase margins of multi-rate sampled-data feedback systems", *International Journal of Control* 44(3), 1986.
- P. Seiler, A. Packard and P. Gahinet, "An introduction to disk margins", *IEEE Control Systems Magazine* 40(5), 2020.
