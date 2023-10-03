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
    x = [Num(Symbolics.variable(name; T = FnType{Tuple{Any},Real}))(t) for name in x]
    u = [Num(Symbolics.variable(name; T = FnType{Tuple{Any},Real}))(t) for name in u]
    y = [Num(Symbolics.variable(name; T = FnType{Tuple{Any},Real}))(t) for name in y]
    u = map(u) do u
        ModelingToolkit.setmetadata(u, ModelingToolkit.VariableInput, true)
    end
    y = map(y) do y
        ModelingToolkit.setmetadata(y, ModelingToolkit.VariableOutput, true)
    end

    Blocks.StateSpace(ssdata(sys)...; x = x0, name)
end


"""
    sconnect(input::Function, sys::T; name)
    sconnect(input::Num,      sys::T; name)

Connect a function `input(t)` to `sys.input`

# Examples:
```julia
sconnect(sin, sys)   # Connect a funciton, assumed to be a function of time
sconnect(sin(t), sys) # Connect a Num
```
"""
function sconnect(
    input::Union{Function, Num},
    sys::T;
    name = Symbol("$(sys.name) with input"),
) where {T<:ModelingToolkit.AbstractTimeDependentSystem}
    @named output = Blocks.RealOutput()
    T(
        [
            sys.input.u ~ (input isa Num ? input : input(t))
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
    outputs;
    kwargs...
)
    named_ss(sys, inputs, outputs; kwargs...).sys # just discard the names
end


"""
    RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractTimeDependentSystem, inputs, outputs; kwargs...)

Convert an `ODESystem` to a `NamedStateSpace` using linearization. `inputs, outputs` are vectors of variables determining the inputs and outputs respectively. See docstring of `ModelingToolkit.linearize` for more info on `kwargs`.

This method automatically converts systems that MTK has failed to produce a proper form for into a proper linear statespace system. Learn more about how that is done here:
https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/#Internals:-Transformation-of-non-proper-models-to-proper-statespace-form

See also [`ModelingToolkit.linearize`](@ref) which is the lower-level function called internally. The functions [`get_named_sensitivity`](@ref), [`get_named_comp_sensitivity`](@ref), [`get_named_looptransfer`](@ref) similarily provide convenient ways to compute sensitivity functions while retaining signal names in the same way as `named_ss`. The corresponding lower-level functions `get_sensitivity`, `get_comp_sensitivity` and `get_looptransfer` are available in ModelingToolkitStandardLibrary.Blocks and are documented in [MTKstdlib: Linear analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/).
"""
function RobustAndOptimalControl.named_ss(
    sys::ModelingToolkit.AbstractTimeDependentSystem,
    inputs,
    outputs;
    kwargs...,
)

    if isa(inputs,  Symbol)
        nu = 1
    else
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
        nu = length(inputs)
    end
    if isa(outputs,  Symbol)
        ny = 1
    else
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
        ny = length(outputs)
    end
    matrices, ssys = ModelingToolkit.linearize(sys, inputs, outputs; kwargs...)
    symstr(x) = Symbol(string(x))
    unames = symstr.(inputs)
    fm(x) = convert(Matrix{Float64}, x)
    if nu > 0 && size(matrices.B, 2) == 2nu
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
        Ce = [fm(matrices.C) zeros(ny, ndu)]
        De = fm(matrices.D[:, 1:nu])
        dsys = dss(Ae, E, Be, Ce, De)
        lsys = ss(RobustAndOptimalControl.DescriptorSystems.dss2ss(dsys)[1])
        # unames = [unames; Symbol.("der_" .* string.(unames))]
        # sys = ss(matrices...)
    else
        lsys = ss(matrices...)
    end
    named_ss(
        lsys;
        x = symstr.(states(ssys)),
        u = unames,
        y = symstr.(outputs),
        name = string(Base.nameof(sys)),
    )
end

for f in [:sensitivity, :comp_sensitivity, :looptransfer]
    fnn = Symbol("get_named_$f")
    fn = Symbol("get_$f")
    @eval function $(fnn)(args...; kwargs...)
        named_sensitivity_function(Blocks.$(fn), args...; kwargs...)
    end
end


"""
    get_named_sensitivity(sys, ap::AnalysisPoint; kwargs...)
    get_named_sensitivity(sys, ap_name::Symbol; kwargs...)

Call [`ModelingToolkitStandardLibrary.Blocks.get_sensitivity`](@ref) while retaining signal names. Returns a `NamedStateSpace` object (similar to [`named_ss`](@ref)).
"""
get_named_sensitivity

"""
    get_named_comp_sensitivity(sys, ap::AnalysisPoint; kwargs...)
    get_named_comp_sensitivity(sys, ap_name::Symbol; kwargs...)

Call [`ModelingToolkitStandardLibrary.Blocks.get_comp_sensitivity`](@ref) while retaining signal names. Returns a `NamedStateSpace` object (similar to [`named_ss`](@ref)).
"""
get_named_comp_sensitivity

"""
    get_named_looptransfer(sys, ap::AnalysisPoint; kwargs...)
    get_named_looptransfer(sys, ap_name::Symbol; kwargs...)

Call [`ModelingToolkitStandardLibrary.Blocks.get_looptransfer`](@ref) while retaining signal names. Returns a `NamedStateSpace` object (similar to [`named_ss`](@ref)).
"""
get_named_looptransfer

function named_sensitivity_function(
    fun,
    sys::ModelingToolkit.AbstractTimeDependentSystem,
    inputs, args...;
    kwargs...,
)

    if isa(inputs,  Symbol)
        nu = 1
    else
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
        nu = length(inputs)
    end
    matrices, ssys = fun(sys, inputs, args...; kwargs...)
    symstr(x) = Symbol(string(x))
    unames = symstr.(inputs)
    fm(x) = convert(Matrix{Float64}, x)
    if nu > 0 && size(matrices.B, 2) == 2nu
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
        Ce = [fm(matrices.C) zeros(ny, ndu)]
        De = fm(matrices.D[:, 1:nu])
        dsys = dss(Ae, E, Be, Ce, De)
        lsys = ss(RobustAndOptimalControl.DescriptorSystems.dss2ss(dsys)[1])
        # unames = [unames; Symbol.("der_" .* string.(unames))]
        # sys = ss(matrices...)
    else
        lsys = ss(matrices...)
    end
    named_ss(
        lsys;
        x = symstr.(states(ssys)),
        u = unames,
        y = unames, #Symbol.("out_" .* string.(inputs)),
        name = string(Base.nameof(sys)),
    )
end

if isdefined(ModelingToolkit, :get_disturbance_system)
    function ModelingToolkit.get_disturbance_system(dist::ModelingToolkit.DisturbanceModel{<:LTISystem})
        ControlSystemsBase.issiso(dist.model) || error("Disturbance model must be SISO")
        Blocks.StateSpace(ssdata(ss(dist.model))..., name=dist.name)
    end
end

"""
    build_quadratic_cost_matrix(linear_sys, ssys::ODESystem, costs::Vector{Pair})

For a system that has been linearized, assemble a quadratic cost matrix (for LQR or Kalman filtering) that penalizes states or outputs of simplified system `ssys` according to the vector of pairs `costs`.

The motivation for this function is that ModelingToolkit does not guarantee
- Which states are selected as states after simplification.
- The order of the states.

The second problem above, the ordering of the states, can be worked around using `reorder_states`, but the first problem cannot be solved by trivial reordering. This function thus accepts an array of costs for a user-selected state realization, and assembles the correct cost matrix for the state realization selected by MTK. To do this, the funciton needs the linearization (`linear_sys`) as well as the simplified system, both of which are outputs of [`linearize`](@ref).

# Arguments:
- `linear_sys`: Output of [`linearize`](@ref), an object containing a property called `C`. This can be a [`ControlSystemsBase.StateSpace`](@ref) or a `NamedTuple` with a field `C`.
- `ssys`: Output of [`linearize`](@ref).
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

Assemble a quadratic cost matrix (for LQR or Kalman filtering) that penalizes states or outputs of system `sys` according to the vector of pairs `costs`.

The motivation for this function is that ModelingToolkit does not guarantee
- Which states are selected as states after simplification.
- The order of the states.

The second problem above, the ordering of the states, can be worked around using `reorder_states`, but the first problem cannot be solved by trivial reordering. This function thus accepts an array of costs for a user-selected state realization, and assembles the correct cost matrix for the state realization selected by MTK. To do this, the funciton performs a linearization between inputs and the cost outputs. The linearization is used to determine the matrix entries belonging to states that are not part of the realization chosen by MTK.

# Arguments:
- `sys`: The system to be linearized (not simplified).
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


function batch_linearize(sys, inputs, outputs, ops::AbstractVector{<:AbstractDict}; t = 0.0,
        allow_input_derivatives = false,
        kwargs...)
    lin_fun, ssys = linearization_function(sys, inputs, outputs; kwargs...)
    lins = map(ops) do op
        linearize(ssys, lin_fun; op, t, allow_input_derivatives)
    end
    lins, ssys
end

"""
    batch_ss(sys, inputs, outputs, ops::AbstractVector{<:AbstractDict};
                t = 0.0,
                allow_input_derivatives = false,
                kwargs...)

Linearize `sys` in multiple operating points `ops::Vector{Dict}`. Returns a vector of `StateSpace` objects and the simplified system.

# Example:
```
using ControlSystemsMTK, ModelingToolkit, RobustAndOptimalControl
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

bounds = getbounds(duffing, states(duffing))
sample_within_bounds((l, u)) = (u - l) * rand() + l
# Create a vector of operating points
ops = map(1:N) do i
    op = Dict(x => sample_within_bounds(bounds[x]) for x in keys(bounds) if isfinite(bounds[x][1]))
end

Ps, ssys = batch_ss(duffing, [u], [y], ops)
w = exp10.(LinRange(-2, 2, 200))
bodeplot(Ps, w)
P = RobustAndOptimalControl.ss2particles(Ps) # convert to a single StateSpace system with `Particles` as coefficients.
bodeplot(P, w) # Should look similar to the one above
```

Let's also do some tuning for the linearized models above
```
function batch_tune(f, Ps)
    f.(Ps)
end

Cs = batch_tune(Ps) do P
    # C, kp, ki, fig, CF = loopshapingPI(P, 6; phasemargin=45)
    C, kp, ki, kd, fig, CF = loopshapingPID(P, 6; Mt=1.3, Tf = 1/100)
    ss(CF)
end

P = RobustAndOptimalControl.ss2particles(Ps)
C = RobustAndOptimalControl.ss2particles(Cs)

nyquistplot(P * C,
            w,
            ylims = (-10, 3),
            xlims = (-5, 10),
            points = true,
            Ms_circles = [1.5, 2],
            Mt_circles = [1.5, 2])

# Fit circles that encircle the Nyquist curve for each frequency
centers, radii = fit_complex_perturbations(P * C, w; relative = false, nominal = :center)
nyquistcircles!(w, centers, radii, ylims = (-4, 1), xlims = (-3, 4))
```

See also [`trajectory_ss`](@ref) and [`fuzz`](@ref).
"""
function batch_ss(args...; kwargs...)
    lins, ssys = batch_linearize(args...; kwargs...)
    [ss(l...) for l in lins], ssys
end

"""
    linsystems, ssys = trajectory_ss(sys, inputs, outputs, sol; t = _max_100(sol.t), fuzzer=nothing, verbose = true, kwargs...)

Linearize `sys` around the trajectory `sol` at times `t`. Returns a vector of `StateSpace` objects and the simplified system.

# Arguments:
- `inputs`: A vector of variables or analysis points.
- `outputs`: A vector of variables or analysis points.
- `sol`: An ODE solution object. This solution must contain the states of the simplified system, accessible through the `idxs` argument like `sol(t, idxs=x)`.
- `t`: Time points along the solution trajectory at which to linearize. The returned array of `StateSpace` objects will be of the same length as `t`.
- `fuzzer`: A function that takes an operating point dictionary and returns an array of "fuzzed" operating points. This is useful for adding noise/uncertainty to the operating points along the trajectory. See [`ControlSystemsMTK.fuzz`](@ref) for such a function.
- `verbose`: If `true`, print warnings for variables that are not found in `sol`.
- `kwargs`: Are sent to the linearization functions.
"""
function trajectory_ss(sys, inputs, outputs, sol; t = _max_100(sol.t), allow_input_derivatives = false, fuzzer = nothing, verbose = true, kwargs...)
    maximum(t) > maximum(sol.t) && @warn("The maximum time in `t`: $(maximum(t)), is larger than the maximum time in `sol.t`: $(maximum(sol.t)).")
    minimum(t) < minimum(sol.t) && @warn("The minimum time in `t`: $(minimum(t)), is smaller than the minimum time in `sol.t`: $(minimum(sol.t)).")
    lin_fun, ssys = linearization_function(sys, inputs, outputs; kwargs...)
    x = states(ssys)
    defs = ModelingToolkit.defaults(sys)
    ops = map(t) do ti
        Dict(x => robust_sol_getindex(sol, ti, x, defs; verbose) for x in x)
    end
    if fuzzer !== nothing
        opsv = map(ops) do op
            fuzzer(op)
        end
        ops = reduce(vcat, opsv)
        t = repeat(t, inner = length(ops) ÷ length(t))
    end
    lins = map(zip(ops, t)) do (op, t)
        linearize(ssys, lin_fun; op, t, allow_input_derivatives)
    end
    (; linsystems = [ss(l...) for l in lins], ssys, ops)
end

"_max_100(t) = length(t) > 100 ? range(extrema(t)..., 100) : t"
_max_100(t) = length(t) > 100 ? range(extrema(t)..., 100) : t

"""
    fuzz(op, p; N = 10, parameters = true, variables = true)

"Fuzz" an operating point `op::Dict` by changing each non-zero value to an uncertain number with multiplicative uncertainty `p`, represented by `N` samples, i.e., `p = 0.1` means that the value is multiplied by a `N` numbers between 0.9 and 1.1.

`parameters` and `variables` indicate whether to fuzz parameters and state variables, respectively.

This function modifies all variables the same way. For more fine-grained control, load the `MonteCarloMeasurements` package and use the `Particles` type directly, followed by `MonteCarloMeasurements.particle_dict2dict_vec(op)`, i.e., the following makes `uncertain_var` uncertain with a 10% uncertainty:
```julia
using MonteCarloMeasurements
op = ModelingToolkit.defaults(sys)
op[uncertain_var] = op[uncertain_var] * Particles(10, Uniform(0.9, 1.1))
ops = MonteCarloMeasurements.particle_dict2dict_vec(op)
batch_ss(model, inputs, outputs, ops)
```
If you have more than one uncertain parameter, it's important to use the same number of particles for all of them (10 in the example above).

To make use of this function in [`trajectory_ss`](@ref), pass something like
```
fuzzer = op -> ControlSystemsMTK.fuzz(op, 0.02; N=10)
```
to fuzz each operating point 10 times with a 2% uncertainty. The resulting number of operating points will increase by 10x.
"""
function fuzz(op, p; N=10, parameters = true, variables = true)
    op = map(collect(keys(op))) do key
        par = ModelingToolkit.isparameter(key)
        val = op[key]
        par && !parameters && return (key => val)
        !par && !variables && return (key => val)
        aval = abs(val)
        uval = issymbolic(val) ? val : iszero(val) ? 0.0 : Particles(N, MonteCarloMeasurements.Uniform(val-p*aval, val+p*aval))
        key => uval
    end |> Dict
    MonteCarloMeasurements.particle_dict2dict_vec(op)
end

MonteCarloMeasurements.vecindex(p::Symbolics.BasicSymbolic,i) = p
issymbolic(x) = x isa Union{Symbolics.Num, Symbolics.BasicSymbolic}

"""
    robust_sol_getindex(sol, ti, x, defs; verbose = true)

Extract symbolic variable `x` from ode solution `sol` at time `ti`. This operation may fail
- If the variable is a dummy derivative that is not present in the solution. In this case, the value is reconstructed by derivative interpolation.
- The var is not present at all, in this case, the default value in `defs` is returned.

# Arguments:
- `sol`: An ODESolution
- `ti`: Time point
- `defs`: A Dict with default values. 
- `verbose`: Print a warning if the variable is not found in the solution.
"""
function robust_sol_getindex(sol, ti, x, defs; verbose = true)
    try
        return sol(ti, idxs=x)
    catch
        n = string((x))
        if occursin("ˍt(", n)
            n = split(n, "ˍt(")[1]
            sp = split(n, '₊')
            varname = sp[end]
            local var
            let t = Symbolics.arguments(Symbolics.unwrap(x))[1]
                @variables var(t)
            end
            ModelingToolkit.@set! var.val.f.name = Symbol(varname)
            namespaces = sp[1:end-1]
            if !isempty(namespaces)
                for ns in reverse(namespaces)
                    var = ModelingToolkit.renamespace(Symbol(ns), var)
                end
            end
            out = sol(ti, Val{1}, idxs=[Num(var)])[]
            verbose && println("Could not find variable $x in solution, returning $(out) obtained through interpolation of $var.")
            return out
        end

        val = get(defs, x, 0.0)
        verbose && println("Could not find variable $x in solution, returning $val.")
        return val
    end
end

maybe_interp(interpolator, x, t) = allequal(x) ? x[1] : interpolator(x, t)

"""
    GainScheduledStateSpace(systems, vt; interpolator, x = zeros((systems[1]).nx), name, u0 = zeros((systems[1]).nu), y0 = zeros((systems[1]).ny))

A linear parameter-varying (LPV) version of [`Blocks.StateSpace`](@ref), implementing the following equations:

```math
\\begin{aligned}
\\dot{x} &= A(v) x + B(v) u \\\\
y        &= C(v) x + D(v) u
\\end{aligned}
```

where `v` is a scalar scheduling variable.

# Arguments:
- `systems`: A vector of `ControlSystemsBase.StateSpace` objects
- `vt`: A vector of breakpoint values for the scheduling variable `v`, this has the same length as `systems`.
- `interpolator`: A constructor `i = interpolator(values, breakpoints)` and returns an interpolator object that can be called like `i(v)` to get the interpolated value at `v`. `LinearInterpolation` from DataInterpolations.jl is a good choice, but a lookup table can also be used.
"""
function GainScheduledStateSpace(systems, vt; interpolator, x = zeros(systems[1].nx), name, u0 = zeros(systems[1].nu), y0 = zeros(systems[1].ny))

    s1 = first(systems)
    (; nx, nu, ny) = s1

    Aint = [maybe_interp(interpolator, getindex.(getfield.(systems, :A), i, j), vt) for i = 1:nx, j = 1:nx]
    Bint = [maybe_interp(interpolator, getindex.(getfield.(systems, :B), i, j), vt) for i = 1:nx, j = 1:nu]
    Cint = [maybe_interp(interpolator, getindex.(getfield.(systems, :C), i, j), vt) for i = 1:ny, j = 1:nx]
    Dint = [maybe_interp(interpolator, getindex.(getfield.(systems, :D), i, j), vt) for i = 1:ny, j = 1:nu]

    @named input = Blocks.RealInput(nin = nu)
    @named scheduling_input = Blocks.RealInput()
    @named output = Blocks.RealOutput(nout = ny)
    @variables x(t)[1:nx]=x [
        description = "State variables of gain-scheduled statespace system $name",
    ]
    @variables v(t) = 0 [
        description = "Scheduling variable of gain-scheduled statespace system $name",
    ]
    
    @variables A(v)[1:nx, 1:nx] = systems[1].A
    @variables B(v)[1:nx, 1:nu] = systems[1].B
    @variables C(v)[1:ny, 1:nx] = systems[1].C
    @variables D(v)[1:ny, 1:nu] = systems[1].D
    A,B,C,D = collect.((A,B,C,D))

    eqs = [
        v ~ scheduling_input.u;
        [A[i] ~ (Aint[i] isa Number ? Aint[i] : Aint[i](v)) for i in eachindex(A)];
        [B[i] ~ (Bint[i] isa Number ? Bint[i] : Bint[i](v)) for i in eachindex(B)];
        [C[i] ~ (Cint[i] isa Number ? Cint[i] : Cint[i](v)) for i in eachindex(C)];
        [D[i] ~ (Dint[i] isa Number ? Dint[i] : Dint[i](v)) for i in eachindex(D)];
        [Differential(t)(x[i]) ~ sum(A[i, k] * x[k] for k in 1:nx) +
                                 sum(B[i, j] * (input.u[j] - u0[j]) for j in 1:nu)
         for i in 1:nx];
        output.u .~ C * x .+ D * (input.u .- u0) .+ y0
    ]
    compose(ODESystem(eqs, t, name = name), [input, output, scheduling_input])
end

"LPVStateSpace is equivalent to GainScheduledStateSpace, see the docs for GainScheduledStateSpace."
const LPVStateSpace = GainScheduledStateSpace


# struct InterpolatorGain{I} <: Function
#     interpolator::I
# end

# (ig::InterpolatorGain)(v) = ig.interpolator(v)

#=
The HammersteinWienerSystem idea will not really pan outsince there is no way of getting the InterpolatorGain be a function of v
We can however use the PartitionedStateSpace to prove stability of the closed-llop gain-scheduled system by moving the maximum gain of the interpolator to either the bottom row or the rightmost column, and then calculating the structured singular value with real structured uncertainty (infinitely time varying)
=#
# function GainScheduledLFT(systems, vt; interpolator)
#     ig(array) = InterpolatorGain(interpolator(array, vt))
#     s1 = first(systems)
#     (; nx, nu, ny) = s1

#     A,B,C,D = ssdata(s1)
#     Aconst = [allequal(getindex.(getfield.(systems, :A), i, j)) for i = 1:nx, j = 1:nx]
#     Bconst = [allequal(getindex.(getfield.(systems, :B), i, j)) for i = 1:nx, j = 1:nu]
#     Cconst = [allequal(getindex.(getfield.(systems, :C), i, j)) for i = 1:ny, j = 1:nx]
#     Dconst = [allequal(getindex.(getfield.(systems, :D), i, j)) for i = 1:ny, j = 1:nu]
#     A = A .* Aconst
#     B1 = B .* Bconst
#     C1 = C .* Cconst
#     D11 = D .* Dconst

#     Mlpv = .! [Aconst Bconst; Cconst Dconst]
#     total_num_lpv = count(!, Aconst) + count(!, Bconst) + count(!, Cconst) + count(!, Dconst)

#     B2 = zeros(nx, total_num_lpv)
#     C2 = zeros(total_num_lpv, nx)
#     D12 = zeros(ny, total_num_lpv)
#     D21 = zeros(total_num_lpv, nu)
#     D22 = zeros(total_num_lpv, total_num_lpv)

#     c = 1
#     interp_gains = Function[]
#     for i = axes(Mlpv, 1), j = axes(Mlpv, 2)
#         # Incoming LPV are B2, D12, D22
#         # Outgoing LPV are C2, D21, D22
#         # LPV for A has B2 as incoming and C2 as outgoing
#         # LPV for B has B2 as incoming and D21 as outgoing
#         # LPV for C has D12 as incoming and C2 as outgoing
#         # LPV for D has D22 as both incoming and outgoing
#         if Mlpv[i, j]
#             if i <= nx # A or B
#                 if j <= nx # A
#                     B2[i, c] = 1
#                     C2[c, j] = 1
#                     push!(interp_gains, ig([s.A[i,j] for s in systems]))
#                 else # B
#                     B2[c, j - nx] = 1
#                     D21[i - ny, c] = 1
#                     push!(interp_gains, ig([s.B[i,j-nx] for s in systems]))
#                 end
#             else # C or D
#                 if j <= nx # C
#                     D12[i-nx, c] = 1
#                     C2[c, j] = 1
#                     push!(interp_gains, ig([s.C[i-nx,j] for s in systems]))
#                 else # D
#                     D22[c, c] = 1
#                     push!(interp_gains, ig([s.D[i-nx,j-nx] for s in systems]))
#                 end
#             end
#             c += 1
#         end
#     end

#     P = ControlSystemsBase.PartitionedStateSpace(A, B1, B2, C1, C2, D11, D12, D21, D22, s1.timeevol)
#     # ControlSystemsBase.HammersteinWienerSystem(P, interp_gains)
# end


"""
    Symbolics.build_function(sys::AbstractStateSpace, args; kwargs)

Build a function that takes parameters and returns a [`StateSpace`](@ref) object (technically a `HeteroStateSpace`).

# Arguments:
- `sys`: A statespace system with (typically) symbolic coefficients.
- `args` and `kwargs`: are passed to the internal call to `build_function` from the Symbolics.jl package.
"""
function Symbolics.build_function(sys::AbstractStateSpace, args...; kwargs...)
    ControlSystemsBase.numeric_type(sys) <: Num || error("Expected a system with symbolic coefficients. Call linearize_symbolic to obtain symbolic jacobians")
    Afun, _ = Symbolics.build_function(sys.A, args...; kwargs...)
    Bfun, _ = Symbolics.build_function(sys.B, args...; kwargs...)
    Cfun, _ = Symbolics.build_function(sys.C, args...; kwargs...)
    Dfun, _ = Symbolics.build_function(sys.D, args...; kwargs...)
    (args...) -> HeteroStateSpace(Afun(args...), Bfun(args...), Cfun(args...), Dfun(args...), sys.timeevol)
end