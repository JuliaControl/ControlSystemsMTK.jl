using ModelingToolkit, OrdinaryDiffEq
@parameters t
D = Differential(t)

"""
    ESC(; k, tau, a, w, wh = 10w, name)

PI-ESC: Proportional-integral extremum seeking controller. 

Assume a system to be controlled on the form
```math
\\dot x = f(x) + g(x)u \\\\
y = h(x)
```
where ``h(x)`` is a cost function to be minimized, the value of which can be observed (the function may in practice be unknown). The PI-ESC finds the steady-state control signal ``u^*`` that minimizes ``y`` by applying a dither signal ``a \\sin(ωt)`` to estimate the derivative of the cost function w.r.t. the input.

# Arguments:
- `k`: Proportional gain
- `tau`: Integral time constant
- `a`: Dither amplitude
- `w`: Dither frequency
- `wh`: High-pass filter frequency, typically 10x larger than `w`.

Ref: "Proportional-integral extremum-seeking control" M. Guay
"""
function ESC(; k, tau, a, w, wh=10*w, name)
    wh > w || error("wh should be (typically 10x) larger than w")
    @parameters k=k tau=tau a=a w=w wh=wh
    @variables v(t)=0 uh(t)=0 u(t)=0 y(t)=0
    eqs = [
        D(v)  ~ -wh*v + y
        D(uh) ~ -1/tau * (-wh^2*v + wh*y)*sin(w*t)
        u     ~ -k/a   * (-wh^2*v + wh*y)*sin(w*t) + uh + a*sin(w*t)
    ]
    ODESystem(eqs, t; name)
end
