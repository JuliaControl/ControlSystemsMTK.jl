using ModelingToolkit: AnalysisPoint
const AP = Union{Symbol, AnalysisPoint}
import ModelingToolkitStandardLibrary.Blocks as Blocks
conn = ModelingToolkit.connect
t = Blocks.t
ModelingToolkit.System(sys::LTISystem; kwargs...) = System(ss(sys); kwargs...)

"""
    ModelingToolkit.System(sys::AbstractStateSpace; name::Symbol, x0 = zeros(sys.nx), x_names, u_names, y_names)

Create an System from `sys::StateSpace`. 

# Arguments:
- `sys`: An instance of `StateSpace` or `NamedStateSpace`.
- `name`: A symbol giving the system a unique name.
- `x0`: Initial state
The arguments below are automatically set if the system is a `NamedStateSpace`.
- `x_names`: A vector of symbols with state names. 
- `u_names`: A vector of symbols with input names. 
- `y_names`: A vector of symbols with output names. 
"""
function ModelingToolkit.System(
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
    uc = [Blocks.RealInput(; name = Symbol(u)) for u in u]
    yc = [Blocks.RealOutput(; name = Symbol(y)) for y in y]
    @named ssblock = Blocks.StateSpace(ssdata(sys)...; x = x0)
    @unpack input, output = ssblock
    systems = [uc; yc; input; output]
    eqs = [
        [uc[i].u ~ input.u[i] for i in 1:length(uc)];
        [yc[i].u ~ output.u[i] for i in 1:length(yc)];
    ]
    extend(System(eqs, t; name, systems), ssblock)
end


numeric(x::Num) = x.val

function ControlSystemsBase.ss(
    sys::ModelingToolkit.AbstractSystem,
    inputs,
    outputs;
    kwargs...
)
    named_ss(sys, inputs, outputs; kwargs...).sys # just discard the names
end

symstr(x) = Symbol(x isa AnalysisPoint ? x.name : string(x))


"""
    RobustAndOptimalControl.named_ss(sys::ModelingToolkit.AbstractSystem, inputs, outputs; descriptor=true, simple_infeigs=true, kwargs...)

Convert an `System` to a `NamedStateSpace` using linearization. `inputs, outputs` are vectors of variables determining the inputs and outputs respectively. See docstring of `ModelingToolkit.linearize` for more info on `kwargs`.

If `descriptor = true` (default), this method automatically converts systems that MTK has failed to produce a proper form for into a proper linear statespace system using the method described here:
https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/#Internals:-Transformation-of-non-proper-models-to-proper-statespace-form
If `descriptor = false`, the system is instead converted to a statespace realization using `sys[:,uinds] + sys[:,duinds]*tf('s')`, which tends to result in a larger realization on which the user may want to call `minreal(sys, tol)` with a carefully selected tolerance.


See also [`ModelingToolkit.linearize`](@ref) which is the lower-level function called internally. The functions [`get_named_sensitivity`](@ref), [`get_named_comp_sensitivity`](@ref), [`get_named_looptransfer`](@ref) similarily provide convenient ways to compute sensitivity functions while retaining signal names in the same way as `named_ss`. The corresponding lower-level functions `get_sensitivity`, `get_comp_sensitivity` and `get_looptransfer` are available in ModelingToolkitStandardLibrary.Blocks and are documented in [MTKstdlib: Linear analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/).
"""
function RobustAndOptimalControl.named_ss(
    sys::ModelingToolkit.AbstractSystem,
    inputs,
    outputs;
    descriptor = true,
    simple_infeigs = true,
    kwargs...,
)

    inputs = vcat(inputs)
    outputs = vcat(outputs)

    inputs = map(inputs) do inp
        if inp isa System
            @variables u(t)
            if u ∈ Set(unknowns(inp))
                inp.u
            else
                error("Input $(inp.name) is an System and not a variable")
            end
        else
            inp
        end
    end
    nu = length(inputs)

    outputs = map(outputs) do out
        if out isa System
            @variables u(t)
            if u ∈ Set(unknowns(out))
                out.u
            else
                error("Outut $(out.name) is an System and not a variable")
            end
        else
            out
        end
    end
    matrices, ssys, xpt = ModelingToolkit.linearize(sys, inputs, outputs; kwargs...)
    unames = symstr.(inputs)
    if nu > 0 && size(matrices.B, 2) == 2nu
        # This indicates that input derivatives are present
        duinds = findall(any(!iszero, eachcol(matrices.B[:, nu+1:end]))) .+ nu
        u2du = (1:nu) .=> duinds # This maps inputs to their derivatives
        lsys = causal_simplification(matrices, u2du; descriptor, simple_infeigs)
    else
        lsys = ss(matrices...)
    end
    pind = [ModelingToolkit.parameter_index(ssys, i) for i in ModelingToolkit.inputs(ssys)]
    x0 = xpt.x
    u0 = [xpt.p[pi] for pi in pind] 
    xu = (; x = x0, u = u0)
    extra = Dict(:operating_point => xu)
    # If simple_infeigs=false, the system might have been reduced and the state names might not match the original system.
    x_names = simple_infeigs ? symstr.(unknowns(ssys)) : [Symbol(string(nameof(sys))*"_x$i") for i in 1:lsys.nx]
    named_ss(
        lsys;
        x = x_names,
        u = unames,
        y = symstr.(outputs),
        name = string(Base.nameof(sys)),
        extra,
    )
end

"""
    causal_simplification(sys, u2duinds::Vector{Pair{Int, Int}}; descriptor=true, simple_infeigs = true)
"""
function causal_simplification(sys, u2duinds::Vector{Pair{Int, Int}}; descriptor=true, simple_infeigs = true)
    fm(x) = convert(Matrix{Float64}, x)
    nx = size(sys.A, 1)
    ny = size(sys.C, 1)
    ndu = length(u2duinds)
    nu = size(sys.B, 2) - ndu
    u_with_du_inds = first.(u2duinds)
    duinds = last.(u2duinds)
    B = sys.B[:, 1:nu]
    B̄ = sys.B[:, duinds]
    D = sys.D[:, 1:nu]
    D̄ = sys.D[:, duinds]
    iszero(fm(D̄)) || error("Nonzero feedthrough matrix from input derivative not supported")
    if descriptor
        Iu = u_with_du_inds .== (1:nu)'
        E = [I(nx) -B̄; zeros(ndu, nx+ndu)]

        Ae = cat(sys.A, -I(ndu), dims=(1,2))
        # Ae[1:nx, nx+1:end] .= B
        Be = [B; Iu]
        Ce = [fm(sys.C) zeros(ny, ndu)]
        De = fm(D)
        dsys = dss(Ae, E, Be, Ce, De)
        return ss(RobustAndOptimalControl.DescriptorSystems.dss2ss(dsys; simple_infeigs)[1])
    else
        ss(sys.A, B, sys.C, D) + ss(sys.A, B̄, sys.C, D̄)*tf('s')
    end
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
    sys::ModelingToolkit.AbstractSystem,
    inputs, args...;
    descriptor = true,
    simple_infeigs = true,
    kwargs...,
)

    inputs = vcat(inputs)
    inputs = map(inputs) do inp
    if inp isa System
            @variables u(t)
            if u ∈ Set(unknowns(inp))
                inp.u
            else
                error("Input $(inp.name) is an System and not a variable")
            end
        else
            inp
        end
    end
    nu = length(inputs)
    matrices, ssys = fun(sys, inputs, args...; kwargs...)
    symstr(x) = Symbol(x isa AnalysisPoint ? x.name : string(x))
    unames = symstr.(inputs)
    fm(x) = convert(Matrix{Float64}, x)
    if nu > 0 && size(matrices.B, 2) == 2nu
        # This indicates that input derivatives are present
        duinds = findall(any(!iszero, eachcol(matrices.B[:, nu+1:end]))) .+ nu
        u2du = (1:nu) .=> duinds # This maps inputs to their derivatives
        lsys = causal_simplification(matrices, u2du; descriptor, simple_infeigs)
    else
        lsys = ss(matrices...)
    end
    x_names = simple_infeigs ? symstr.(unknowns(ssys)) : [Symbol(string(nameof(sys))*"_x$i") for i in 1:lsys.nx]
    named_ss(
        lsys;
        x = x_names,
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
    build_quadratic_cost_matrix(linear_sys, ssys::System, costs::Vector{Pair})

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
function build_quadratic_cost_matrix(matrices::NamedTuple, ssys::System, costs::AbstractVector{<:Pair})
    x = ModelingToolkit.unknowns(ssys)
    y = ModelingToolkit.outputs(ssys)
    # y = getproperty.(ModelingToolkit.observed(ssys), :lhs)
    nx = length(x)
    new_Cs = map(costs) do (xi, ci)
        i = findfirst(isequal(xi), x)
        if i !== nothing
            sqrt(ci) .* ((1:nx)' .== i)
        else # not a state, get output instead
            i = findfirst(isequal(xi), y)
            i === nothing && error("$xi is neither a state variable nor an output of the system")
            sqrt(ci) .* matrices.C[i, :]
        end
    end
    C = reduce(vcat, new_Cs)
    C'C
end

"""
    build_quadratic_cost_matrix(sys::System, inputs::Vector, costs::Vector{Pair}; kwargs...)

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
function build_quadratic_cost_matrix(sys::System, inputs::AbstractVector, costs::AbstractVector{<:Pair}; kwargs...)
    matrices, ssys, extras = ModelingToolkit.linearize(sys, inputs, first.(costs); kwargs...)
    x = ModelingToolkit.unknowns(ssys)
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
    lin_fun, ssys = linearization_function(sys, inputs, outputs; op=ops[1], kwargs...)
    lins_ops = map(ops) do op
        linearize(ssys, lin_fun; op, t, allow_input_derivatives)
    end
    lins = first.(lins_ops)
    resolved_ops = last.(lins_ops)
    lins, ssys, resolved_ops
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
@variables y(t)

D = Differential(t)

eqs = [D(x) ~ v
       D(v) ~ -k * x - k3 * x^3 - c * v + 10u
       y ~ x]


@named duffing = System(eqs, t)

bounds = getbounds(duffing, unknowns(duffing))
sample_within_bounds((l, u)) = (u - l) * rand() + l
# Create a vector of operating points
N = 10
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
    lins, ssys, resolved_ops = batch_linearize(args...; kwargs...)
    named_linsystems = map(lins) do l
        # Convert to a NamedStateSpace with the same names as the original system
        named_ss(ss(l.A, l.B, l.C, l.D); name = string(Base.nameof(ssys)), x = symstr.(unknowns(ssys)))
    end
    named_linsystems, ssys, resolved_ops
end

# function unnamespace(ap)
#     map(ap.outputs) do out
#         ap_name = ModelingToolkit.SymbolicIndexingInterface.getname(out.u) 
#         new_name = join(ModelingToolkit.namespace_hierarchy(ap_name)[2:end], Symbolics.NAMESPACE_SEPARATOR)
#         Symbolics.rename(ap.input.u, Symbol(new_name))
#     end
# end

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
    # NOTE: we call linearization_funciton twice :( The first call is to get x=unknowns(ssys), the second call provides the operating points.
    # lin_fun, ssys = linearization_function(sys, inputs, outputs; warn_initialize_determined = false, kwargs...)
    lin_fun, ssys = linearization_function(sys, inputs, outputs; warn_empty_op = false, warn_initialize_determined = false, kwargs...)
    x = unknowns(ssys)

    # TODO: The value of the output (or input) of the input analysis points should be mapped to the perturbation vars
    perturbation_vars = ModelingToolkit.inputs(ssys)
    # original_inputs = reduce(vcat, unnamespace(ap) for ap in vcat(inputs)) # assuming all inputs are analysis points for now

    input_names = reduce(vcat, getproperty.(ap.outputs, :u) for ap in vcat(inputs)) 
    output_names = reduce(vcat, ap.input.u for ap in vcat(outputs)) 

    op_nothing = Dict(unknowns(sys) .=> nothing) # Remove all defaults present in the original system
    defs = ModelingToolkit.defaults(sys)
    ops = map(t) do ti
        opsol = Dict(x => robust_sol_getindex(sol, ti, x, defs; verbose) for x in x)
        # opsolu = Dict(new_u => robust_sol_getindex(sol, ti, u, defs; verbose) for (new_u, u) in zip(perturbation_vars, original_inputs))
        merge(op_nothing, opsol)
    end
    if fuzzer !== nothing
        opsv = map(ops) do op
            fuzzer(op)
        end
        ops = reduce(vcat, opsv)
        t = repeat(t, inner = length(ops) ÷ length(t))
    end
    lin_fun, ssys = linearization_function(sys, inputs, outputs; op=ops[1], kwargs...)#, initialization_abstol=1e-1, initialization_reltol=1e-1, kwargs...) # initializealg=ModelingToolkit.SciMLBase.NoInit()
    # Main.lin_fun = lin_fun
    # Main.op1 = ops[1]
    # Main.ops = ops 
    # equations(lin_fun.prob.f.initialization_data.initializeprob.f.sys)
    # observed(lin_fun.prob.f.initialization_data.initializeprob.f.sys)
    lins_ops = map(zip(ops, t)) do (op, t)
        linearize(ssys, lin_fun; op, t, allow_input_derivatives)
        # linearize(sys, inputs, outputs; op, t, allow_input_derivatives) # useful for debugging
    end
    lins = first.(lins_ops)
    resolved_ops = last.(lins_ops)
    named_linsystems = map(lins) do l
        # Convert to a NamedStateSpace with the same names as the original system
        named_ss(ss(l.A, l.B, l.C, l.D); name = string(Base.nameof(sys)), x = symstr.(unknowns(ssys)), u = symstr.(input_names), y = symstr.(output_names))
    end
    (; linsystems = named_linsystems, ssys, ops, resolved_ops)
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

See example usage in the [gain-scheduling example](https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/batch_linearization/#Gain-scheduling).

# Arguments:
- `systems`: A vector of `ControlSystemsBase.StateSpace` objects
- `vt`: A vector of breakpoint values for the scheduling variable `v`, this has the same length as `systems`.
- `interpolator`: A constructor `i = interpolator(values, breakpoints)` and returns an interpolator object that can be called like `i(v)` to get the interpolated value at `v`. `LinearInterpolation` from DataInterpolations.jl is a good choice, but a lookup table can also be used.

# Connectors
- `input` of type `RealInput` connects to ``u``.
- `output` of type `RealOutput` connects to ``y``.
- `scheduling_input` of type `RealInput` connects to ``v``.
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
    @variables v(t) [
        description = "Scheduling variable of gain-scheduled statespace system $name",
    ]
    
    @variables A(t)[1:nx, 1:nx] = systems[1].A
    @variables B(t)[1:nx, 1:nu] = systems[1].B
    @variables C(t)[1:ny, 1:nx] = systems[1].C
    @variables D(t)[1:ny, 1:nu] = systems[1].D
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
        collect(output.u .~ C * x .+ D * (input.u .- u0) .+ y0)
    ]
    compose(System(eqs, t, name = name), [input, output, scheduling_input])
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
