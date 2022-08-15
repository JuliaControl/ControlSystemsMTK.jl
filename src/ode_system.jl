import ModelingToolkitStandardLibrary.Blocks as Blocks
conn = ModelingToolkit.connect
t = Blocks.t
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
function ModelingToolkit.ODESystem(sys::AbstractStateSpace;
    name::Symbol,
    x0 = zeros(sys.nx),
    x_names = [Symbol("x$i") for i in 1:sys.nx],
    u_names = sys.nu == 1 ? [:u] : [Symbol("u$i") for i in 1:sys.nu],
    y_names = sys.ny == 1 ? [:y] : [Symbol("y$i") for i in 1:sys.ny],
)
    ControlSystems.isdiscrete(sys) && error("Discrete systems not yet supported due to https://github.com/SciML/ModelingToolkit.jl/issues?q=is%3Aopen+is%3Aissue+label%3Adiscrete-time")
    A,B,C,D = ssdata(sys)
    nx,ny,nu = sys.nx, sys.ny, sys.nu
    x = [Num(Symbolics.variable(name; T=FnType{Tuple{Any},Real}))(t) for name in x_names]
    u = [Num(Symbolics.variable(name; T=FnType{Tuple{Any},Real}))(t) for name in u_names]
    y = [Num(Symbolics.variable(name; T=FnType{Tuple{Any},Real}))(t) for name in y_names]
    u = map(u) do u
        ModelingToolkit.setmetadata(u, ModelingToolkit.VariableInput, true)
    end
    y = map(y) do y
        ModelingToolkit.setmetadata(y, ModelingToolkit.VariableOutput, true)
    end

    osys = Blocks.StateSpace(ssdata(sys)...; x_start = x0, name)
end

function ModelingToolkit.ODESystem(sys::NamedStateSpace;
    name::Symbol,
    kwargs...
)
    @unpack x_names, u_names, y_names = sys
    ODESystem(sys.sys; x_names, u_names, y_names, name, kwargs...)
end


"""
    sconnect(input, sys::T; name)
"""
function sconnect(input, sys::T; name=Symbol("$(sys.name) with input")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    T([
        conn(input.output, sys.input)
    ], t; systems=[sys, input], name)
end

"""
    sconnect(input::Function, sys::T; name)

Connect a function `input(t)` to `sys.input`
"""
function sconnect(input::Function, sys::T; name=Symbol("$(sys.name) with input")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    @named output = Blocks.RealOutput()
    T([
        sys.input.u ~ input(t)
        output.u ~ sys.output.u
    ], t; systems=[sys, output], name)
end

"""
    sconnect(sys1::T, sys2::T; name)

Connect systems in series, equivalent to `sys2*sys1` or `series(sys1, sys2)` in ControlSystems.jl terminology
"""
function sconnect(sys1::T, sys2::T; name=Symbol("$(sys1.name)*$(sys2.name)")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    @named output = Blocks.RealOutput() # TODO: missing size
    @named input = Blocks.RealInput() # TODO: missing size
    T([
        conn(input, sys2.input)
        conn(output, sys1.output)
        conn(sys2.output, sys1.input)
    ], t; name, systems=[sys1, sys2, output, input])
end

"""
    G = ControlSystems.feedback(loopgain::T; name)

Form the feedback-interconnection
\$G = L/(1+L)\$

The system `G` will be a new system with `input` and `output` connectors.
"""
function ControlSystems.feedback(loopgain::T; name=Symbol("feedback $(loopgain.name)")) where T <: ModelingToolkit.AbstractTimeDependentSystem
    add = Blocks.Add(k1=1, k2=-1, name=:feedback)
    @named input = Blocks.RealInput()
    @named output = Blocks.RealOutput()
    T([
        input.u ~ add.input1.u
        output.u ~ loopgain.output.u
        conn(loopgain.output, add.input2)
        conn(add.output, loopgain.input)
    ], t; systems=[input, output, loopgain, add], name)
end

function Base.:(*)(s1::T, s2::T) where T <: ModelingToolkit.AbstractTimeDependentSystem
    name = Symbol(string(s1.name)*"_"*string(s2.name))
    @named input = Blocks.RealInput()
    @named output = Blocks.RealOutput()
    eqs = [
        conn(s1.input, s2.output)
        output.u ~ s1.output.u
    ]
    systems=[output, s1, s2]
    if any(s.name == :input for s in s2.systems)
        push!(eqs, input.u ~ s2.input.u)
        push!(systems, input)
    end
    T(eqs, t; systems, name)
end


numeric(x::Num) = x.val


function ControlSystems.ss(sys::ModelingToolkit.AbstractTimeDependentSystem, inputs, outputs)
    named_ss(sys, inputs, outputs).sys # just discard the names
end

inputs(sys) = filter(s->ModelingToolkit.isinput(s), states(sys))
outputs(sys) = filter(s->ModelingToolkit.isoutput(s), states(sys))

"""
    RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractTimeDependentSystem, inputs, outputs; kwargs...)

Convert an `ODESystem` to a `NamedStateSpace`. `inputs, outputs` are vectors of variables determining the inputs and outputs respectively. See docstring of `ModelingToolkit.linearize` for more info on `kwargs`, reproduced below.

$(@doc(ModelingToolkit.linearize))
"""
function RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractTimeDependentSystem, inputs, outputs; kwargs...)

    inputs = map(inputs) do inp
        if inp isa ODESystem
            @variables u(t)
            if u ∈ Set(states(inp))
                inp.u
            else
                error("Input $(inp.name) is an ODESystem and not a variable")
            end
        else
            inp
        end
    end
    outputs = map(outputs) do out
        if out isa ODESystem
            @variables u(t)
            if u ∈ Set(states(out))
                out.u
            else
                error("Outut $(out.name) is an ODESystem and not a variable")
            end
        else
            out
        end
    end
    matrices, ssys = ModelingToolkit.linearize(sys, inputs, outputs; kwargs...)
    symstr(x) = Symbol(string(x))
    named_ss(ss(matrices...); x=symstr.(states(ssys)), u=symstr.(inputs), y=symstr.(outputs))
end
