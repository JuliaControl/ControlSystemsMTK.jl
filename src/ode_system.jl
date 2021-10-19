ModelingToolkit.ODESystem(sys::LTISystem; kwargs...) = ODESystem(ss(sys); kwargs...)

"""
    ModelingToolkit.ODESystem(sys::AbstractStateSpace{Continuous}; name::Symbol, x0 = zeros(sys.nx), x_names, u_names, y_names)

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
    x = [Num(Variable{FnType{Tuple{Any},Real}}(name))(t) for name in x_names]
    u = [Num(Variable{FnType{Tuple{Any},Real}}(name))(t) for name in u_names]
    y = [Num(Variable{FnType{Tuple{Any},Real}}(name))(t) for name in y_names]
    Dₜ = Differential(t)
    eqs = [
        Dₜ.(x) .~ A*x .+ B*u
        y      .~ C*x .+ D*u
    ]
    ODESystem(eqs; name, defaults = Dict(x .=> x0))
end

function ModelingToolkit.ODESystem(sys::NamedStateSpace{Continuous};
    name::Symbol,
    x0 = zeros(sys.nx)
)
    @unpack x_names, u_names, y_names = sys
    ODESystem(sys.sys; x_names, u_names, y_names)
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


function ControlSystems.ss(sys::ODESystem, inputs, outputs)
    named_ss(sys, inputs, outputs).sys # just discard the names
end

# This approach is very brittle 
"""
    RobustAndOptimalControl.named_ss(sys::ODESystem, u, y)

Convert an `ODESystem` to a `NamedStateSpace`. `u,y` are vectors of variables determining the inputs and outputs respectively. If the input `u` contains a `@register`ed function, apply the function to the corresponding symbolic argument, e.g.,
```
@parameter t
@register input(t)
named_ss(sys, input(t), y)
```
"""
# function RobustAndOptimalControl.named_ss(sys::ODESystem, u, y)
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

function RobustAndOptimalControl.named_ss(sys::ODESystem)
    # u isa AbstractVector || (u = [u])
    # y isa AbstractVector || (y = [y])
    u = controls(sys)
    y = observed(sys)
    x = states(sys)
    
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
    named_ss(ss(A,B,C,D); x=ff.(x), u=ff.(u), y=ff.(ylhs))
end
