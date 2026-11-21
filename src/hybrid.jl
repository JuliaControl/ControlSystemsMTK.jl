# Linearization of hybrid models, that is, models in which clocked (discrete-time) partitions
# are connected to continuous-time equations.

# A discrete compilation pass, such as `SynchToolkit.compile_lustre`, among the
# `additional_passes` indicates that the model has clocked partitions. With such a pass,
# `ModelingToolkit.linearize` succeeds but its result omits the clocked partitions.
function has_discrete_compile_pass(kwargs)
    passes = get(kwargs, :additional_passes, nothing)
    passes === nothing && return false
    return any(ModelingToolkit.discrete_compile_pass, passes)
end

# `linearize_hybrid` compiles each partition itself, the discrete compilation pass is not used
function without_discrete_compile_pass(kwargs)
    has_discrete_compile_pass(kwargs) || return kwargs
    passes = filter(!ModelingToolkit.discrete_compile_pass, collect(kwargs[:additional_passes]))
    return (; kwargs..., additional_passes = passes)
end

function check_no_discrete_compile_pass(sys, fname, kwargs)
    has_discrete_compile_pass(kwargs) || return nothing
    throw(ArgumentError("A discrete compilation pass was supplied in `additional_passes`, which indicates that the model `$(nameof(sys))` contains clocked partitions. `$fname` does not support such models: the linearization would omit the clocked partitions. Use `named_ss` or the `get_named_*` functions, which assemble the sampled-data system of a hybrid model, or `ModelingToolkit.linearize_hybrid`."))
end

# Call `linearization()`, the linearization that does not account for clocked partitions.
# Without a discrete compilation pass, ModelingToolkit rejects a model with clocked
# partitions; unless the caller chose `hybrid = false`, `fallback()` is then called.
function linearize_or_fallback(linearization, fallback, hybrid)
    try
        return linearization()
    catch e
        (hybrid === nothing && e isa ModelingToolkit.HybridSystemNotSupportedException) || rethrow()
    end
    return fallback()
end

# Call the hybrid counterpart of the linearization function `fun`
function hybrid_linearization(fun, sys, args...; kwargs...)
    kwargs = without_discrete_compile_pass(kwargs)
    if !isdefined(ModelingToolkit, :linearize_hybrid)
        throw(ArgumentError("The model `$(nameof(sys))` contains clocked partitions. Linearization of such a model requires `ModelingToolkit.linearize_hybrid`, which the installed version of ModelingToolkit does not provide. A linearization that omits the clocked partitions is obtained with `hybrid = false`."))
    end
    if fun === ModelingToolkit.linearize
        ModelingToolkit.linearize_hybrid(sys, args...; DEFAULT_LINEARIZE_KWARGS..., kwargs...)
    else
        fun(sys, args...; hybrid = true, DEFAULT_LINEARIZE_KWARGS..., kwargs...)
    end
end

# The scalar elements of the variables `signals`, with array-valued variables expanded
function scalar_signals(signals)
    out = Any[]
    for s in signals
        s = Symbolics.unwrap(s)
        if Symbolics.symbolic_type(s) == Symbolics.ArraySymbolic()
            append!(out, Symbolics.unwrap.(vec(collect(s))))
        else
            push!(out, s)
        end
    end
    return out
end

# Position of each user signal in the concatenation of the signals of all partitions.
# `user[i]` holds the user signals of partition `i` and `offsets[i]` the number of signals of
# the partitions preceding it. Within one partition the user signals appear in the order
# given by the user; when they are distributed over several partitions, the order is
# recovered by matching against `signals`, which is only possible for variables.
function user_signal_indices(user, offsets, signals, kind)
    occupied = findall(!isempty, user)
    if length(occupied) <= 1
        return reduce(vcat, (offsets[i] .+ eachindex(user[i]) for i in occupied); init = Int[])
    end
    if signals === nothing || any(s -> s isa Union{Symbol, AnalysisPoint}, signals)
        throw(ArgumentError("The $kind of the linearization belong to several clock partitions and are specified through analysis points, in which case their order cannot be recovered. Specify the $kind as variables, or use `ModelingToolkit.linearize_hybrid` and assemble the partitions manually."))
    end
    return map(scalar_signals(signals)) do s
        for i in occupied
            k = findfirst(isequal(s), user[i])
            k === nothing || return offsets[i] + k
        end
        error("The signal $s was not found among the $kind of the hybrid linearization.")
    end
end

"""
    hybrid_ss(hl; descriptor = true, simple_infeigs = true, balance = false, big = false)

Assemble the partitions of the `ModelingToolkit.HybridLinearization` `hl` into a single
`StateSpace` system from all user-specified inputs to all user-specified outputs. Returns
`(sys, uinds, yinds, x_names)` where `sys` has the inputs and outputs of all partitions,
including the signals crossing clock boundaries, and `x_names` is `nothing` if the state
realization of a partition was changed during the assembly.

If the model has clocked partitions, the result is a discrete-time system: the continuous
partition is discretized with zero-order hold at the sample interval of the clock, and the
connections across the clock boundaries are closed. All clocked partitions are required to
be on the same periodic clock, and the continuous partition must not have direct feedthrough
from a held signal to a sampled signal.
"""
function hybrid_ss(hl; descriptor = true, simple_infeigs = true, balance = false, big = false)
    partitions = hl.partitions
    ci = hl.continuous_index
    discrete = [p for (i, p) in enumerate(partitions) if i != ci]
    if isempty(discrete)
        Ts = nothing
    else
        clocks = unique(p.clock for p in discrete)
        if length(clocks) > 1 || any(p -> p.Ts === nothing, discrete)
            throw(ArgumentError("The assembly of a hybrid linearization into a single system requires all clocked partitions to be on the same periodic clock, but the model has the clocks $(join(string.(clocks), ", ")). Use `ModelingToolkit.linearize_hybrid` and assemble the partitions manually."))
        end
        Ts = discrete[1].Ts
    end
    named_states = true
    systems = map(enumerate(partitions)) do (i, p)
        i == ci || return ss(p.A, p.B, p.C, p.D, Ts)
        nu = length(p.inputs)
        if nu > 0 && size(p.B, 2) == 2nu
            # Input derivatives are present, see `named_ss`
            u2du = [j => j + nu for j in 1:nu]
            P = causal_simplification(p, u2du; descriptor, simple_infeigs, big, balance)
            named_states = false
        else
            P = ss(p.A, p.B, p.C, p.D)
        end
        Ts === nothing ? P : c2d(P, Ts)
    end
    if ci !== nothing
        # A sampled signal with direct feedthrough from a held signal depends on the value
        # that was held before the tick, which the instantaneous connection below does not
        # represent
        held = unique(c.input for c in hl.connections if c.to == ci)
        sampled = unique(c.output for c in hl.connections if c.from == ci)
        if !iszero(systems[ci].D[sampled, held])
            throw(ArgumentError("The continuous-time partition has direct feedthrough from a held signal to a sampled signal. The assembly of such a model into a single system is not supported. Use `ModelingToolkit.linearize_hybrid` and assemble the partitions manually."))
        end
    end
    A = cat((s.A for s in systems)...; dims = (1, 2))
    B = cat((s.B for s in systems)...; dims = (1, 2))
    C = cat((s.C for s in systems)...; dims = (1, 2))
    D = cat((s.D for s in systems)...; dims = (1, 2))
    uoffsets = cumsum([0; [s.nu for s in systems[1:end-1]]])
    yoffsets = cumsum([0; [s.ny for s in systems[1:end-1]]])
    # The connections across clock boundaries are u = K*y + w, where w are the inputs of the
    # assembled system
    K = zeros(size(B, 2), size(C, 1))
    for c in hl.connections
        K[uoffsets[c.to] + c.input, yoffsets[c.from] + c.output] = 1
    end
    E = lu(I - K * D; check = false)
    issuccess(E) || throw(ArgumentError("The connections across the clock boundaries form an algebraic loop, the assembled system is not well posed."))
    KC = E \ (K * C)
    Acl = A + B * KC
    Bcl = B / E
    Ccl = C + D * KC
    Dcl = D / E
    sys = Ts === nothing ? ss(Acl, Bcl, Ccl, Dcl) : ss(Acl, Bcl, Ccl, Dcl, Ts)
    x_names = named_states ? reduce(vcat, (symstr.(p.unknowns) for p in partitions); init = Symbol[]) : nothing
    if x_names !== nothing && !(allunique(x_names) && length(x_names) == sys.nx)
        x_names = nothing
    end
    return sys, uoffsets, yoffsets, x_names
end

function named_ss_hybrid(hl, sys, inputs, outputs, unames, ynames; kwargs...)
    lsys, uoffsets, yoffsets, x_names = hybrid_ss(hl; kwargs...)
    partitions = hl.partitions
    uinds = user_signal_indices([p.inputs[1:p.nu_user] for p in partitions], uoffsets, inputs, "inputs")
    yinds = user_signal_indices([p.outputs[1:p.ny_user] for p in partitions], yoffsets, outputs, "outputs")
    lsys = lsys[yinds, uinds]
    if x_names === nothing
        x_names = [Symbol(string(nameof(sys))*"_x$i") for i in 1:lsys.nx]
    end
    nsys = named_ss(
        lsys;
        x = x_names,
        u = unames,
        y = ynames,
        name = string(Base.nameof(sys)),
    )
    RobustAndOptimalControl.set_extra!(nsys, :hybrid_linearization, hl)
    nsys
end
