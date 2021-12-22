ModelingToolkit.ODESystem(sys::LTISystem; kwargs...) = ODESystem(ss(sys); kwargs...)

"""
    ModelingToolkit.ODESystem(sys::AbstractStateSpace; name::Symbol, x0 = zeros(sys.nx), x_names, u_names, y_names)

Create an ODESystem from `sys::StateSpace`. 

# Arguments:
- `sys`: An instance of `StateSpace` or `NamedStateSpace`.
- `name`: A symbol giving the system a unique name.
- `x0`: Initial state
The arguments below are automatically set if the system is a `NamedStateSpace`.
- `x_names`: A vector of symbols with state names. 
- `u_names`: A vector of symbols with input names. 
- `y_names`: A vector of symbols with output names. 
"""
function system_creator(constructor, sys::AbstractStateSpace;
    name::Symbol,
    x0 = zeros(sys.nx),
    x_names = [Symbol("x$i") for i in 1:sys.nx],
    u_names = sys.nu == 1 ? [:u] : [Symbol("u$i") for i in 1:sys.nu],
    y_names = sys.ny == 1 ? [:y] : [Symbol("y$i") for i in 1:sys.ny],
)
    A,B,C,D = ssdata(sys)
    nx,ny,nu = sys.nx, sys.ny, sys.nu
    @parameters t
    x = [Num(Symbolics.variable(name; T=FnType{Tuple{Any},Real}))(t) for name in x_names]
    u = [Num(Symbolics.variable(name; T=FnType{Tuple{Any},Real}))(t) for name in u_names] # TODO: should be input=true
    y = [Num(Symbolics.variable(name; T=FnType{Tuple{Any},Real}))(t) for name in y_names] # TODO: should be output=true
    # @show typeof(u)
    # u = map(u) do u
    #     ModelingToolkit.toinput(u.val).f
    # end
    # y = map(y) do y
    #     ModelingToolkit.tooutput(y.val).f
    # end
    Dₜ = if ControlSystems.isdiscrete(sys)
        # A = A - I # Due to difference operator instead of time shift operator https://github.com/SciML/ModelingToolkit.jl/issues/1307
        # NOTE: The difference operator appears to be a time-shift operator
        Difference(t; dt=sys.Ts)
    else
        Differential(t)
    end
    eqs = [
        Dₜ.(x) .~ A*x .+ B*u
        y      .~ C*x .+ D*u
        ]
    constructor(eqs, t; name, defaults = Dict([x .=> x0; y .=> C*x0; u .=> 0])) # NOTE: the initial values for y and u are required since `structural_simplify` fails on discrete systems
end

function ModelingToolkit.ODESystem(sys::AbstractStateSpace, args...; kwargs...)
    system_creator(ODESystem, sys, args...; kwargs...)
end

function ModelingToolkit.DiscreteSystem(sys::AbstractStateSpace{<:Discrete}, args...; kwargs...)
    system_creator(DiscreteSystem, sys, args...; kwargs...)
end

function ModelingToolkit.ODESystem(sys::NamedStateSpace{Continuous};
    name::Symbol,
    kwargs...
)
    @unpack x_names, u_names, y_names = sys
    ODESystem(sys.sys; x_names, u_names, y_names, name, kwargs...)
end

function ModelingToolkit.ODESystem(sys::NamedStateSpace{<:Discrete};
    name::Symbol,
    kwargs...
)
    @unpack x_names, u_names, y_names = sys
    ODESystem(sys.sys; x_names, u_names, y_names, name, kwargs...)
end

function sconnect(input, sys::T; name=Symbol("$(sys.name) with input")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    @parameters t
    @variables u(t) y(t)
    T([
        u ~ input
        sys.u ~ input
        y ~ sys.y
    ], t; systems=[sys], name)
end

function sconnect(input::Function, sys::T; name=Symbol("$(sys.name) with input")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    @parameters t
    @variables u(t) y(t)
    T([
        sys.u ~ input(u)
        y ~ sys.y
    ], t; systems=[sys], name)
end

"connect output of one sys to input of other"
function sconnect(sys1::T, sys2::T; name=Symbol("$(sys1.name)*$(sys2.name)")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    @parameters t
    @variables u(t) y(t)
    T([
        u ~ sys1.u
        sys1.y ~ sys2.u
        y ~ sys2.y
    ], t; systems=[sys1, sys2], name)
end

"form feedback interconnection, i.e., input is `r-y`"
function ControlSystems.feedback(loopgain::T, ref; name=Symbol("feedback $(loopgain.name)")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    @parameters t
    @variables u(t) y(t)
    T([
        u ~ ref
        ref - loopgain.y ~ loopgain.u
        y ~ loopgain.y
    ], t; systems=[loopgain], name)
end

function ControlSystems.feedback(loopgain::T, ref::T; name=Symbol("feedback $(loopgain.name)*$(ref.name)")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    @parameters t
    @variables u(t) y(t)
    T([
        u ~ ref.u
        ref.y - loopgain.y ~ loopgain.u
        y ~ loopgain.y
    ], t; systems=[loopgain, ref], name)
end

numeric(x::Num) = x.val


function ControlSystems.ss(sys::ModelingToolkit.AbstractTimeDependentSystem, inputs, outputs)
    named_ss(sys, inputs, outputs).sys # just discard the names
end

inputs(sys) = filter(s->ModelingToolkit.isinput(s), states(sys))
outputs(sys) = filter(s->ModelingToolkit.isoutput(s), states(sys))

"""
    RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractTimeDependentSystem, u, y)

Convert an `ODESystem` to a `NamedStateSpace`. `u,y` are vectors of variables determining the inputs and outputs respectively. If the input `u` contains a `@register`ed function, apply the function to the corresponding symbolic argument, e.g.,
```
@parameter t
@register input(t)
named_ss(sys, input(t), y)
```
"""
function RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractTimeDependentSystem; Ts = nothing, numeric=false,
    u = controls(sys),
    y = observed(sys),
    x = states(sys)
)

    isempty(u) && error("No control inputs specified.")
    isempty(y) && error("No outputs specified.")
    
    equation_order = map(equations(sys)) do eq # dynamics equations and states are not in the same order

    end
    A = calculate_jacobian(sys)
    B = calculate_control_jacobian(sys) 
    yrhs = [e.rhs for e in y]
    ylhs = [e.lhs for e in y]
    C = jacobian(yrhs, x)
    D = jacobian(yrhs, u)
    # D = 0
    @assert size(C, 1) == length(y) "C matrix of wrong size: $C, aeqs: $aeqs"
    @assert size(D) == (length(y), length(u)) "D matrix of wrong size: $D"
    symstr(x) = Symbol(string(x))
    ff(x) = symstr(x.f)
    if numeric
        sym2val = ModelingToolkit.defaults(sys)
        A = substitute.(A, Ref(sym2val)) .|> ControlSystemsMTK.numeric .|> Float64
        B = substitute.(B, Ref(sym2val)) .|> ControlSystemsMTK.numeric .|> Float64
        C = substitute.(C, Ref(sym2val)) .|> ControlSystemsMTK.numeric .|> Float64
        D = substitute.(D, Ref(sym2val)) .|> ControlSystemsMTK.numeric .|> Float64
    end
    timeevol = if Ts === nothing
        sys isa ODESystem || error("Sample time Ts must be provided for discrete-time systems")
        Continuous()
    else
        Discrete(Ts)
    end
    named_ss(ss(A,B,C,D,timeevol); x=ff.(x), u=ff.(u), y=ff.(ylhs))
end

# This approach is very brittle 
# function RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractTimeDependentSystem, u, y)
#     u isa AbstractVector || (u = [u])
#     y isa AbstractVector || (y = [y])
#     eqs     = equations(sys)
#     # find all differential terms and filter x based on those
#     lhss    = getproperty.(eqs, :lhs)
#     dx      = filter(isdifferential, lhss)
#     x       = [dx.arguments[] for dx in dx]
#     diffeqs = [e.rhs for e in eqs if Symbolics.is_derivative(e.lhs)]
#     sy      = Set(y)
#     # aeqsl    = [e.rhs for e in eqs if e.lhs ∈ sy] # This fails if the equation is not on the right form. Below is a failed attempt at solving this.
#     # aeqsr    = [e.lhs for e in eqs if e.rhs ∈ sy]
#     # aeqs     = [aeqsl; aeqsr]
#     # might have to go through all observed variables and if an y is in there, make substitution
#     o = observed(sys)
#     olhs = Set([o.lhs for o in o])
#     # y = map(y) do y # go through all outputs and replace if they are observed
#     #     if y isa Equation # this can happen if y is observed
#     #         y = y.lhs
#     #     end
#     #     if y ∈ olhs
#     #         i = findfirst(isequal(y), olhs
#     #         return o[i].rhs
#     #     end
#     #     y
#     # end
#     eqvars = map(eqs) do eq
#         s1 = Set(Symbolics.get_variables(eq.lhs))
#         s2 = Set(Symbolics.get_variables(eq.rhs))
#         [s1; s2]
#     end
#     # start by going through observed and get those y equations, if a particular y is not observed, fall back to other method. 
#     # the aeqs are 1) observed, 2) solvefor y
#     aeqs = map(y) do y
#         y isa Equation && (return y.rhs) # in this case y is already on the required form. This happens if y ∈ observed(sys)
#         for (i,eqv) in enumerate(eqvars)
#             y ∈ eqv || continue
#             return solve_for(eqs[i], y)
#         end
#     end
#     @show aeqs
#     @show jacobian(aeqs, states(sys)) # This line fails, probably due to namespacing issues 
#     A = jacobian(diffeqs, x) .|> numeric
#     B = jacobian(diffeqs, u) .|> numeric
#     C = jacobian(aeqs, x)    .|> numeric
#     D = jacobian(aeqs, u)    .|> numeric
#     @assert size(C, 1) == length(y) "C matrix of wrong size: $C, aeqs: $aeqs"
#     @assert size(D) == (length(y), length(u)) "D matrix of wrong size: $D"
#     symstr(x) = Symbol(string(x))
#     named_ss(ss(A,B,C,D); x_names=symstr.(x), u_names=symstr.(u), y_names=symstr.(y))
# end