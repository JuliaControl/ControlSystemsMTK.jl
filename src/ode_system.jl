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
function ModelingToolkit.ODESystem(
    sys::AbstractStateSpace;
    name::Symbol,
    x0 = zeros(sys.nx),
    x = ControlSystemsBase.state_names(sys),
    u = ControlSystemsBase.input_names(sys),
    y = ControlSystemsBase.output_names(sys),
)
    ControlSystemsBase.isdiscrete(sys) && error(
        "Discrete systems not yet supported due to https://github.com/SciML/ModelingToolkit.jl/issues?q=is%3Aopen+is%3Aissue+label%3Adiscrete-time",
    )
    A, B, C, D = ssdata(sys)
    nx, ny, nu = sys.nx, sys.ny, sys.nu
    x = [Num(Symbolics.variable(name; T = FnType{Tuple{Any},Real}))(t) for name in x]
    u = [Num(Symbolics.variable(name; T = FnType{Tuple{Any},Real}))(t) for name in u]
    y = [Num(Symbolics.variable(name; T = FnType{Tuple{Any},Real}))(t) for name in y]
    u = map(u) do u
        ModelingToolkit.setmetadata(u, ModelingToolkit.VariableInput, true)
    end
    y = map(y) do y
        ModelingToolkit.setmetadata(y, ModelingToolkit.VariableOutput, true)
    end

    osys = Blocks.StateSpace(ssdata(sys)...; x_start = x0, name)
end

"""
    sconnect(input, sys::T; name)
"""
function sconnect(
    input,
    sys::T;
    name = Symbol("$(sys.name) with input"),
) where {T<:ModelingToolkit.AbstractTimeDependentSystem}
    T([conn(input.output, sys.input)], t; systems = [sys, input], name)
end

"""
    sconnect(input::Function, sys::T; name)

Connect a function `input(t)` to `sys.input`
"""
function sconnect(
    input::Function,
    sys::T;
    name = Symbol("$(sys.name) with input"),
) where {T<:ModelingToolkit.AbstractTimeDependentSystem}
    @named output = Blocks.RealOutput()
    T(
        [
            sys.input.u ~ input(t)
            output.u ~ sys.output.u
        ],
        t;
        systems = [sys, output],
        name,
    )
end

"""
    sconnect(sys1::T, sys2::T; name)

Connect systems in series, equivalent to `sys2*sys1` or `series(sys1, sys2)` in ControlSystems.jl terminology
"""
function sconnect(
    sys1::T,
    sys2::T;
    name = Symbol("$(sys1.name)*$(sys2.name)"),
) where {T<:ModelingToolkit.AbstractTimeDependentSystem}
    @named output = Blocks.RealOutput() # TODO: missing size
    @named input = Blocks.RealInput() # TODO: missing size
    T(
        [
            conn(input, sys2.input)
            conn(output, sys1.output)
            conn(sys2.output, sys1.input)
        ],
        t;
        name,
        systems = [sys1, sys2, output, input],
    )
end

"""
    G = ControlSystemsBase.feedback(loopgain::T; name)

Form the feedback-interconnection
\$G = L/(1+L)\$

The system `G` will be a new system with `input` and `output` connectors.
"""
function ControlSystemsBase.feedback(
    loopgain::T;
    name = Symbol("feedback $(loopgain.name)"),
) where {T<:ModelingToolkit.AbstractTimeDependentSystem}
    add = Blocks.Add(k1 = 1, k2 = -1, name = :feedback)
    @named input = Blocks.RealInput()
    @named output = Blocks.RealOutput()
    T(
        [
            input.u ~ add.input1.u
            output.u ~ loopgain.output.u
            conn(loopgain.output, add.input2)
            conn(add.output, loopgain.input)
        ],
        t;
        systems = [input, output, loopgain, add],
        name,
    )
end

function Base.:(*)(s1::T, s2::T) where {T<:ModelingToolkit.AbstractTimeDependentSystem}
    name = Symbol(string(s1.name) * "_" * string(s2.name))
    @named input = Blocks.RealInput()
    @named output = Blocks.RealOutput()
    eqs = [
        conn(s1.input, s2.output)
        output.u ~ s1.output.u
    ]
    systems = [output, s1, s2]
    if any(s.name == :input for s in s2.systems)
        push!(eqs, input.u ~ s2.input.u)
        push!(systems, input)
    end
    T(eqs, t; systems, name)
end


numeric(x::Num) = x.val


function ControlSystemsBase.ss(
    sys::ModelingToolkit.AbstractTimeDependentSystem,
    inputs,
    outputs,
)
    named_ss(sys, inputs, outputs).sys # just discard the names
end

inputs(sys) = filter(s -> ModelingToolkit.isinput(s), states(sys))
outputs(sys) = filter(s -> ModelingToolkit.isoutput(s), states(sys))

"""
    RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractTimeDependentSystem, inputs, outputs; kwargs...)

Convert an `ODESystem` to a `NamedStateSpace`. `inputs, outputs` are vectors of variables determining the inputs and outputs respectively. See docstring of `ModelingToolkit.linearize` for more info on `kwargs`, reproduced below.

$(@doc(ModelingToolkit.linearize))
"""
function RobustAndOptimalControl.named_ss(
    sys::ModelingToolkit.AbstractTimeDependentSystem,
    inputs,
    outputs;
    kwargs...,
)

    if inputs isa Symbol
        outputs isa Symbol || throw(ArgumentError("inputs and outputs must be either both symbols or both vectors of symbols"))
        nu = ny = 1
    else # map symbols to symbolic variables
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
        nu = length(inputs)
        ny = length(outputs)
    end
    matrices, ssys = ModelingToolkit.linearize(sys, inputs, outputs; kwargs...)
    symstr(x) = Symbol(string(x))
    unames = symstr.(inputs)
    if size(matrices.B, 2) == 2nu
        nx = size(matrices.A, 1)
         # This indicates that input derivatives are present
        duinds = findall(any(!iszero, eachcol(matrices.B[:, nu+1:end])))
        B̄ = matrices.B[:, duinds .+ nu]
        ndu = length(duinds)
        B = matrices.B[:, 1:nu]
        Iu = duinds .== (1:nu)'
        E = [I(nx) -B̄; zeros(ndu, nx+ndu)]

        Ae = cat(matrices.A, -I(ndu), dims=(1,2))
        Be = [B; Iu]
        Ce = [matrices.C zeros(ny, ndu)]
        De = matrices.D[:, 1:nu]
        dsys = dss(Ae, E, Be, Ce, De)
        sys = ss(RobustAndOptimalControl.DescriptorSystems.dss2ss(dsys)[1])
        # unames = [unames; Symbol.("der_" .* string.(unames))]
        # sys = ss(matrices...)
    else
        sys = ss(matrices...)
    end
    named_ss(
        sys;
        x = symstr.(states(ssys)),
        u = unames,
        y = symstr.(outputs),
    )
end

if isdefined(ModelingToolkit, :get_disturbance_system)
    function ModelingToolkit.get_disturbance_system(dist::ModelingToolkit.DisturbanceModel{<:LTISystem})
        ControlSystemsBase.issiso(dist.model) || error("Disturbance model must be SISO")
        Blocks.StateSpace(ssdata(ss(dist.model))..., name=dist.name)
    end
end

"""
    build_quadratic_cost_matrix(matrices::NamedTuple, ssys::ODESystem, costs::Vector{Pair})

For a system that has been linearized, assemble a quadratic cost matrix (for LQR or Kalman filtering) that penalizes states or outputs of simplified system `ssys` accoring to the vector of pairs `costs`.

The motivation for this function is that ModelingToolkit does not guarantee
- Which states are selected as states after simplification.
- The order of the states.

The second problem above, the ordering of the states, can be worked around using `reorder_states`, but the first problem cannot be solved by trivial reordering. This function thus accepts an array of costs for a user-selected state realization, and assembles the correct cost matrix for the state realization selected by MTK. To do this, the funciton needs the linearization (`matrices`) as well as the simplified system, both of which are outputs of `linearize`.

# Arguments:
- `matrices`: Output of `linearize`, an object containing a property called `C`.
- `ssys`: Output of `linearize`.
- `costs`: A vector of pairs
"""
function build_quadratic_cost_matrix(matrices::NamedTuple, ssys::ODESystem, costs::AbstractVector{<:Pair})
    x = ModelingToolkit.states(ssys)
    y = ModelingToolkit.outputs(ssys)
    nx = length(x)
    new_Cs = map(costs) do (xi, ci)
        i = findfirst(isequal(xi), x)
        if i !== nothing
            sqrt(ci) .* ((1:nx)' .== i)
        else # not a state, get output instead
            i = findfirst(isequal(xi), y)
            i === nothing && error("$xi is neither a state nor an output")
            sqrt(ci) .* matrices.C[i, :]
        end
    end
    C = reduce(vcat, new_Cs)
    C'C
end

"""
    build_quadratic_cost_matrix(sys::ODESystem, inputs::Vector, costs::Vector{Pair}; kwargs...)

Assemble a quadratic cost matrix (for LQR or Kalman filtering) that penalizes states or outputs of system `sys` accoring to the vector of pairs `costs`.

The motivation for this function is that ModelingToolkit does not guarantee
- Which states are selected as states after simplification.
- The order of the states.

The second problem above, the ordering of the states, can be worked around using `reorder_states`, but the first problem cannot be solved by trivial reordering. This function thus accepts an array of costs for a user-selected state realization, and assembles the correct cost matrix for the state realization selected by MTK. To do this, the funciton needs the linearization (`matrices`) as well as the simplified system, both of which are outputs of `linearize`.

# Arguments:
- `inputs`: A vector of variables that are to be considered controlled inputs for the LQR controller.
- `costs`: A vector of pairs.
"""
function build_quadratic_cost_matrix(sys::ODESystem, inputs::AbstractVector, costs::AbstractVector{<:Pair}; kwargs...)
    matrices, ssys = ModelingToolkit.linearize(sys, inputs, first.(costs); kwargs...)
    x = ModelingToolkit.states(ssys)
    y = ModelingToolkit.outputs(ssys)
    nx = length(x)
    new_Cs = map(costs) do (xi, ci)
        i = findfirst(isequal(xi), x)
        if i !== nothing
            sqrt(ci) .* ((1:nx)' .== i)
        else # not a state, get output instead
            i = findfirst(isequal(xi), y)
            i === nothing && error("$xi is neither a state nor an output")
            sqrt(ci) .* matrices.C[i, :]
        end
    end
    C = reduce(vcat, new_Cs)
    C'C
end