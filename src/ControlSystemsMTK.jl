module ControlSystemsMTK

#=
Ideas: All connections handled by ModelingToolkit.
Names: 
- handled either by named system, or directly in constructor to ODESystem. 
Functions: 
- Give me linear system from [u1, u3] to [qm, a]
If the linearization of a full system produces a named system, one could implement getindex for vectors of names and obtain the desired transfer functions.


Another idea: use modelingtoolkitize/build_function to obtain a function that can be differentiated with ForwardDiff or FD, https://discourse.julialang.org/t/differentialequations-control-systems-and-linearization/31178/6


A third idea: just use named systems with named indexing to obtain any system you want.

=#
using ModelingToolkit, ControlSystems
using ControlSystems: ssdata, AbstractStateSpace, Continuous, nstates, noutputs, ninputs
import ModelingToolkit: ODESystem
using DifferentialEquations
using ModelingToolkit.Symbolics

export sconnect, feedback, ODESystem

ModelingToolkit.ODESystem(sys::LTISystem; kwargs...) = ODESystem(ss(sys); kwargs...)

"Create an ODESystem from ControlSystems.StateSpace"
function ModelingToolkit.ODESystem(sys::AbstractStateSpace{Continuous};
    name::Symbol,
    x0 = zeros(sys.nx),
    x_names = [Symbol("x$i") for i in 1:sys.nx],
    u_names = sys.nu == 1 ? [:u] : [Symbol("u$i") for i in 1:sys.nu],
    y_names = sys.ny == 1 ? [:y] : [Symbol("y$i") for i in 1:sys.ny],
)
    A,B,C,D = ssdata(sys)
    nx,ny,nu = sys.nx, sys.ny, sys.nu
    # ny == nu == 1 || @warn("MIMO systems are poorly supported for now https://github.com/SciML/ModelingToolkit.jl/issues/605")
    @parameters t
    x = [Num(Variable{ModelingToolkit.FnType{Tuple{Any},Real}}(name))(t) for name in x_names]
    u = [Num(Variable{ModelingToolkit.FnType{Tuple{Any},Real}}(name))(t) for name in u_names]
    y = [Num(Variable{ModelingToolkit.FnType{Tuple{Any},Real}}(name))(t) for name in y_names]
    Dₜ = Differential(t)
    eqs = [
        Dₜ.(x) .~ A*x .+ B*u
        y      .~ C*x .+ D*u
    ]
    ODESystem(eqs; name, defaults = Dict(x .=> x0))
end

"connect input to ODESystem"
function sconnect(input, sys::ODESystem; name=Symbol("$(sys.name) with input"))
    @parameters t
    @variables u(t) y(t)
    ODESystem([
            u ~ input
            sys.u ~ input
            y ~ sys.y
        ], t; systems=[sys], name)
end

function sconnect(input::Function, sys::ODESystem; name=Symbol("$(sys.name) with input"))
    @parameters t
    @variables u(t) y(t)
    ODESystem([
            sys.u ~ input(u)
            y ~ sys.y
        ], t; systems=[sys], name)
end

"connect output of one sys to input of other"
function sconnect(sys1::ODESystem, sys2::ODESystem; name=Symbol("$(sys1.name)*$(sys2.name)"))
    @parameters t
    @variables u(t) y(t)
    ODESystem([
            u ~ sys1.u
            sys1.y ~ sys2.u
            y ~ sys2.y
        ], t; systems=[sys1, sys2], name)
end

"form feedback interconnection, i.e., input is `r-y`"
function ControlSystems.feedback(loopgain::ODESystem, ref; name=Symbol("feedback $(loopgain.name)"))
    @parameters t
    @variables u(t) y(t)
    ODESystem([
            u ~ ref
            ref - loopgain.y ~ loopgain.u
            y ~ loopgain.y
        ], t; systems=[loopgain], name)
end

function ControlSystems.feedback(loopgain::ODESystem, ref; name=Symbol("feedback $(loopgain.name)"))
    @parameters t
    @variables u(t) y(t)
    # @variables y(t)
    ODESystem([
            u ~ ref
            ref - loopgain.y ~ loopgain.u
            y ~ loopgain.y
        ], t; systems=[loopgain], name)
end

function ControlSystems.feedback(loopgain::ODESystem, ref::ODESystem; name=Symbol("feedback $(loopgain.name)*$(ref.name)"))
    @parameters t
    @variables u(t) y(t)
    ODESystem([
            u ~ ref.u
            ref.y - loopgain.y ~ loopgain.u
            y ~ loopgain.y
        ], t; systems=[loopgain, ref], name)
end

numeric(x::Num) = x.val

# This approach is very brittle 
function ControlSystems.ss(sys::ODESystem, inputs, outputs)
    inputs, outputs = Set.((inputs, outputs))
    x       = ModelingToolkit.states(sys) # will contain x, and sometimes u, y
    o       = ModelingToolkit.observed(sys) # equations
    x       = [x; getproperty.(o, :lhs)]
    u       = filter(s->s ∈ inputs, x)
    y       = filter(s->s ∈ outputs, x)
    x       = setdiff(setdiff(x, u), y) # remove states that are not inputs or outputs
    # u = inputs
    eqs     = equations(sys)
    diffeqs = [e.rhs for e in eqs if Symbolics.is_derivative(e.lhs)]
    aeqs    = [e.rhs for e in eqs if e.lhs ∈ outputs]
    A       = Symbolics.jacobian(diffeqs, x) .|> numeric
    B       = Symbolics.jacobian(diffeqs, u) .|> numeric
    C       = Symbolics.jacobian(aeqs, x)    .|> numeric
    D       = Symbolics.jacobian(aeqs, u)    .|> numeric
    ss(A,B,C,D)
end


end
