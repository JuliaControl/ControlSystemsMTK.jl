# Linear analysis of multi-rate hybrid models: periodic model on the event grid, lifting over
# the common period, harmonic frequency functions and margins.

"""
    MultiratePort

A signal of a [`MultirateLinearization`](@ref).

# Fields
- `name::Symbol`: The name of the signal.
- `times::Vector{Float64}`: The instants within one period at which the signal has a value. An input is held between these instants, an output is observed at them.
- `channels::Vector{Int}`: The indices of the corresponding channels of the lifted system, one per instant.
- `partition::Int`: The index of the clock partition of the signal in the `HybridLinearization`.
"""
struct MultiratePort
    name::Symbol
    times::Vector{Float64}
    channels::Vector{Int}
    partition::Int
end

"""
    MultirateLinearization

A linear periodic model of a hybrid system with several clocks, created by
[`multirate_linearization`](@ref).

The model is defined on the event grid: the sorted set of all instants within one common
period `Tc` at which a clock ticks or a signal of the continuous-time partition is updated or
observed. Event `k` maps the state before the event to the state before the next event,

```math
\\begin{aligned}
x_{k+1} &= A_k x_k + B_k u \\\\
y_k &= C_k x_k + D_k u
\\end{aligned}
```

where `u` contains the lifted inputs of the whole period and `y_k` the outputs that are
observed at event `k`. The state consists of the state of the continuous-time partition, the
state variables of the clocked partitions, one hold variable per clocked signal that is read by
another partition, and one hold variable per input of the continuous-time partition.
Each partition is represented by one set of Jacobians, evaluated at the operating point and
the time of the linearization; a model that depends explicitly on time is not represented.

# Fields
- `A, B, C, D`: The matrices of each event.
- `output_channels`: The lifted output channels observed at each event.
- `Tc`: The common period.
- `t0`: The instant at which the period starts.
- `event_times`: The instants of the events, in `[t0, t0 + Tc)`.
- `x`: The names of the state variables.
- `inputs`, `outputs`: The signals, see [`MultiratePort`](@ref).
- `hl`: The `ModelingToolkit.HybridLinearization` the model was created from.

See [`lifted_ss`](@ref), [`multipliers`](@ref), [`fundamental_response`](@ref),
[`harmonic_response`](@ref), [`performance_frequency_gain`](@ref),
[`robust_frequency_gain`](@ref) and [`multirate_margins`](@ref).
"""
struct MultirateLinearization
    A::Vector{Matrix{Float64}}
    B::Vector{Matrix{Float64}}
    C::Vector{Matrix{Float64}}
    D::Vector{Matrix{Float64}}
    output_channels::Vector{Vector{Int}}
    Tc::Float64
    t0::Float64
    event_times::Vector{Float64}
    x::Vector{Symbol}
    inputs::Vector{MultiratePort}
    outputs::Vector{MultiratePort}
    hl::Any
end

function Base.show(io::IO, ::MIME"text/plain", ml::MultirateLinearization)
    println(io, "MultirateLinearization with common period ", ml.Tc, ", ", length(ml.event_times), " events and state dimension ", length(ml.x))
    for (kind, ports) in (("Inputs", ml.inputs), ("Outputs", ml.outputs))
        println(io, kind, ":")
        for p in ports
            println(io, "  ", p.name, ": ", length(p.times), " instant", length(p.times) == 1 ? "" : "s", " per period")
        end
    end
end

# The common period of periods `Ts`: their least common multiple. The ratios of the periods
# are required to be rational with denominators of at most `maxden`.
function common_period(Ts; tol = 1e-9, maxden = 10_000)
    T1 = first(Ts)
    ratios = map(Ts) do T
        r = rationalize(Int, T / T1; tol)
        if denominator(r) > maxden || numerator(r) > maxden
            throw(ArgumentError("The clock periods $(join(string.(Ts), ", ")) are not commensurate: no common period exists. `single_rate_ss` provides an approximation for such models."))
        end
        r
    end
    Tc = T1 * lcm(numerator.(ratios)...) / gcd(denominator.(ratios)...)
    # The ticks of every clock are required to repeat after the common period
    for T in Ts
        abs(round(Tc / T) * T - Tc) <= tol * Tc || throw(ArgumentError("The clock periods $(join(string.(Ts), ", ")) are not commensurate within the tolerance. `single_rate_ss` provides an approximation for such models."))
    end
    Tc
end

# The events of one period. `clocks` maps a partition to the period and phase of its
# instants. The instants are grouped into events independently of `t0`: instants that are
# closer than `atol`, also across the boundary of the period, form one event. Returns the
# times of the events in `[t0, t0 + Tc)` in chronological order, for each event the
# partitions with an instant at the event and the index of that instant, and for each
# partition the times of its instants in chronological order.
function event_schedule(clocks, Tc, t0, atol)
    instants = Tuple{Float64, Int}[]
    for (i, (T, phase)) in clocks
        τ0 = mod(phase, T)
        T - τ0 <= atol && (τ0 = 0.0)
        append!(instants, [(τ0 + m * T, i) for m in 0:round(Int, Tc / T)-1])
    end
    sort!(instants)
    groups = Vector{Tuple{Float64, Int}}[]
    for (τ, i) in instants
        if isempty(groups) || τ - groups[end][1][1] > atol
            push!(groups, Tuple{Float64, Int}[])
        end
        push!(groups[end], (τ, i))
    end
    if length(groups) > 1 && groups[1][1][1] + Tc - groups[end][1][1] <= atol
        # The last instants of the period coincide with the first instants of the next one
        prepend!(groups[1], [(τ - Tc, i) for (τ, i) in pop!(groups)])
    end
    for g in groups
        allunique(last.(g)) || throw(ArgumentError("Two instants of the same clock are closer than the tolerance `tol*Tc`."))
    end
    # Shift every event to the first period that starts at t0
    shifts = [Tc * ceil((t0 - g[end][1] - atol) / Tc) for g in groups]
    order = sortperm([g[end][1] + sh for (g, sh) in zip(groups, shifts)])
    event_times = Float64[]
    event_ticks = Dict{Int, Int}[]
    ticks = Dict(i => Float64[] for i in keys(clocks))
    for q in order
        push!(event_times, groups[q][end][1] + shifts[q])
        push!(event_ticks, Dict{Int, Int}())
        for (τ, i) in groups[q]
            push!(ticks[i], τ + shifts[q])
            event_ticks[end][i] = length(ticks[i])
        end
    end
    event_times, event_ticks, ticks
end

# Compute the outputs of the clocked partitions `ticking` at an event. `known(i, j)` returns
# the row for input `j` of partition `i` if it does not depend on an output of a partition
# that ticks at this event, and `nothing` otherwise. An output that depends directly, through
# the feedthrough matrix, on an output computed at the same event is computed after it.
function clocked_outputs(ticking, systems, source, known, base, width)
    ynew = Dict(i => zeros(systems[i].ny, width) for i in ticking)
    done = Dict(i => falses(systems[i].ny) for i in ticking)
    # The row of input `j` of partition `i`, or `nothing` if it is not yet available
    function input_row(i, j)
        row = known(i, j)
        row === nothing || return row
        from, o = source[(i, j)]
        done[from][o] ? ynew[from][o, :] : nothing
    end
    progress = true
    while progress
        progress = false
        for i in ticking, o in 1:systems[i].ny
            done[i][o] && continue
            row = base(i)[o, :]
            ready = true
            for j in 1:systems[i].nu
                d = systems[i].D[o, j]
                iszero(d) && continue
                u = input_row(i, j)
                u === nothing && (ready = false; break)
                row = row + d * u
            end
            ready || continue
            ynew[i][o, :] = row
            done[i][o] = true
            progress = true
        end
    end
    all(all, values(done)) || throw(ArgumentError("The outputs of the clocked partitions that tick at the same instant depend directly on each other in a cycle; the model is not causal."))
    ynew
end

"""
    multirate_linearization(hl; oversample = 1, t0 = 0.0, continuous_rate = :fastest, inputs = nothing, outputs = nothing, tol = 1e-9)
    multirate_linearization(sys, inputs, outputs; oversample, t0, continuous_rate, tol, kwargs...)

Create a [`MultirateLinearization`](@ref): a linear periodic model of a hybrid system whose
clocked partitions are on periodic clocks with commensurate periods and arbitrary phases.
`hl` is a `ModelingToolkit.HybridLinearization`, the result of
`ModelingToolkit.linearize_hybrid`, or of `get_sensitivity`, `get_comp_sensitivity` and
`get_looptransfer` with `hybrid = true`. The second method calls `linearize_hybrid` with
`kwargs`.

The operations at an event are, in this order:
1. Inputs of the continuous-time partition that are specified by the user take their new
   values, if the event is an instant of their grid.
2. The signals of the continuous-time partition that are sampled are evaluated, with the
   values that the clocked partitions held before the event. A sampled signal with direct
   feedthrough from a held signal thus contains the previously held value, and a sampled
   signal with direct feedthrough from a user input contains its new value.
3. The clocked partitions that tick compute their outputs from their state before the
   event, and then update their state. A partition that reads a clocked signal that is
   computed at the same instant, by another partition or by itself, reads the new value.
   Outputs that depend directly on each other in a cycle are rejected.
4. The held values are updated, and outputs of the continuous-time partition that are
   specified by the user are observed, if the event is an instant of their grid.
5. The continuous-time partition is integrated exactly to the next event with constant
   inputs.

# Keyword arguments
- `continuous_rate`: Inputs and outputs of the continuous-time partition are held and observed on a uniform grid. With `:fastest`, the interval of the grid is the smallest clock period divided by `oversample`; with `:period`, it is the common period.
- `oversample`: An integer that refines the grid of the signals of the continuous-time partition, which makes the behavior between the ticks of the clocks visible. The ticks of the clocks are not changed.
- `t0`: The instant at which the period starts. The characteristic multipliers do not depend on `t0`, the input-output behavior does.
- `inputs`, `outputs`: The variables that were used as inputs and outputs of the linearization. They are required to recover the order of the signals if these belong to several clock partitions.
- `tol`: Relative tolerance for the recognition of commensurate periods and of coinciding instants: instants closer than `tol*Tc` are treated as one event.

The linearization is required to be performed in a stationary operating point of the whole
model. This is not verified.
"""
function multirate_linearization(hl; oversample::Integer = 1, t0 = 0.0, continuous_rate = :fastest,
        inputs = nothing, outputs = nothing, input_names = nothing, output_names = nothing, tol = 1e-9,
        descriptor = true, simple_infeigs = true, balance = descriptor && !simple_infeigs, big = false)
    partitions = hl.partitions
    np = length(partitions)
    ci = something(hl.continuous_index, 0)
    clocked = [i for i in 1:np if i != ci]
    isempty(clocked) && throw(ArgumentError("The model has no clocked partitions."))
    for i in clocked
        partitions[i].Ts === nothing && throw(ArgumentError("The clock $(partitions[i].clock) is not a periodic clock with a fixed sample interval."))
    end
    continuous_rate in (:fastest, :period) || throw(ArgumentError("`continuous_rate` must be `:fastest` or `:period`."))
    oversample >= 1 || throw(ArgumentError("`oversample` must be a positive integer."))
    periods = [partitions[i].Ts for i in clocked]
    Tc = common_period(periods; tol)
    atol = tol * Tc

    # Linear systems of the partitions
    named_states = true
    systems = map(1:np) do i
        p = partitions[i]
        if i == ci
            nu = length(p.inputs)
            if nu > 0 && size(p.B, 2) == 2nu
                # Input derivatives are present, see `named_ss`
                named_states = false
                causal_simplification(p, [j => j + nu for j in 1:nu]; descriptor, simple_infeigs, big, balance)
            else
                ss(p.A, p.B, p.C, p.D)
            end
        else
            ss(p.A, p.B, p.C, p.D, p.Ts)
        end
    end

    # Instants of the clocks and of the signals of the continuous-time partition
    clocks = Dict{Int, Tuple{Float64, Float64}}(i => (partitions[i].Ts, partitions[i].clock.phase) for i in clocked)
    if ci != 0 && (partitions[ci].nu_user > 0 || partitions[ci].ny_user > 0)
        clocks[ci] = (continuous_rate === :period ? Tc : minimum(periods) / oversample, 0.0)
    end
    event_times, event_ticks, ticks = event_schedule(clocks, Tc, t0, atol)
    ci != 0 && !haskey(ticks, ci) && (ticks[ci] = Float64[])
    K = length(event_times)

    # Signals, in the order of the partitions, and their order as specified by the user
    uports = [(i, j) for i in 1:np for j in 1:partitions[i].nu_user]
    yports = [(i, j) for i in 1:np for j in 1:partitions[i].ny_user]
    uoffsets = cumsum([0; [p.nu_user for p in partitions[1:end-1]]])
    yoffsets = cumsum([0; [p.ny_user for p in partitions[1:end-1]]])
    uperm = user_signal_indices([p.inputs[1:p.nu_user] for p in partitions], uoffsets, inputs, "inputs")
    yperm = user_signal_indices([p.outputs[1:p.ny_user] for p in partitions], yoffsets, outputs, "outputs")
    uports, yports = uports[uperm], yports[yperm]
    unames = something(input_names, [symstr(partitions[i].inputs[j]) for (i, j) in uports])
    ynames = something(output_names, [symstr(partitions[i].outputs[j]) for (i, j) in yports])
    function make_ports(ports, names)
        offset = 0
        map(zip(ports, names)) do ((i, _), name)
            channels = collect(offset .+ eachindex(ticks[i]))
            offset += length(channels)
            MultiratePort(Symbol(name), ticks[i], channels, i)
        end
    end
    inports = make_ports(uports, unames)
    outports = make_ports(yports, ynames)
    NU = sum(p -> length(p.channels), inports; init = 0)
    uchannel = Dict((uports[q], r) => inports[q].channels[r] for q in eachindex(uports) for r in eachindex(inports[q].channels))
    ychannel = Dict((yports[q], r) => outports[q].channels[r] for q in eachindex(yports) for r in eachindex(outports[q].channels))

    # State layout
    n = 0
    xidx = map(systems) do s
        idx = n .+ (1:s.nx)
        n += s.nx
        idx
    end
    held = unique((c.from, c.output) for c in hl.connections if c.from != ci)
    hidx = Dict(h => (n += 1) for h in held)
    nuc = ci == 0 ? 0 : partitions[ci].nu_user
    huidx = n .+ (1:nuc)
    n += nuc
    x_names = Symbol[]
    for (i, p) in enumerate(partitions)
        if i == ci && !named_states || length(p.unknowns) != systems[i].nx
            append!(x_names, [Symbol("partition$(i)_x$j") for j in 1:systems[i].nx])
        else
            append!(x_names, symstr.(p.unknowns))
        end
    end
    append!(x_names, [Symbol("hold(", symstr(partitions[i].outputs[o]), ")") for (i, o) in held])
    ci != 0 && append!(x_names, [Symbol("hold(", symstr(partitions[ci].inputs[j]), ")") for j in 1:nuc])
    allunique(x_names) || (x_names = [Symbol("x$i") for i in 1:n])

    # The source of each input of each partition that crosses a clock boundary
    source = Dict((c.to, c.input) => (c.from, c.output) for c in hl.connections)

    A = Matrix{Float64}[]; B = Matrix{Float64}[]; C = Matrix{Float64}[]; D = Matrix{Float64}[]
    output_channels = Vector{Int}[]
    unit(q) = (e = zeros(1, n + NU); e[n + q] = 1; e)
    for k in 1:K
        te = event_times[k]
        # Each row of `R` expresses a state variable in the state before the event and the
        # lifted inputs
        R = [Matrix{Float64}(I, n, n) zeros(n, NU)]
        ticking = [i for i in clocked if haskey(event_ticks[k], i)]
        tick_index(i) = event_ticks[k][i]
        analysis = haskey(event_ticks[k], ci)
        if analysis
            for j in 1:nuc
                R[huidx[j], :] = unit(uchannel[((ci, j), tick_index(ci))])
            end
        end
        # Inputs of the continuous-time partition: user inputs followed by held signals
        function continuous_inputs()
            Uc = zeros(systems[ci].nu, n + NU)
            Uc[1:nuc, :] = R[huidx, :]
            for j in nuc+1:systems[ci].nu
                haskey(source, (ci, j)) && (Uc[j, :] = R[hidx[source[(ci, j)]], :])
            end
            Uc
        end
        yc = ci == 0 ? zeros(0, n + NU) : systems[ci].C * R[xidx[ci], :] + systems[ci].D * continuous_inputs()
        xnew = Dict{Int, Matrix{Float64}}()
        rows = Matrix{Float64}[]
        channels = Int[]
        # The row of input `j` of the clocked partition `i` if it does not depend on an
        # output of a partition that ticks at this event: a user input, a sampled signal of
        # the continuous-time partition, or a held signal of a partition that does not tick
        function known(i, j)
            if j <= partitions[i].nu_user
                vec(unit(uchannel[((i, j), tick_index(i))]))
            elseif !haskey(source, (i, j))
                zeros(n + NU)
            else
                from, o = source[(i, j)]
                from == ci ? yc[o, :] : from in ticking ? nothing : R[hidx[(from, o)], :]
            end
        end
        ynew = clocked_outputs(ticking, systems, source, known, i -> systems[i].C * R[xidx[i], :], n + NU)
        for i in ticking
            Ui = zeros(systems[i].nu, n + NU)
            for j in 1:systems[i].nu
                row = known(i, j)
                Ui[j, :] = row === nothing ? ynew[source[(i, j)][1]][source[(i, j)][2], :] : row
            end
            xnew[i] = systems[i].A * R[xidx[i], :] + systems[i].B * Ui
            for j in 1:partitions[i].ny_user
                push!(rows, ynew[i][j:j, :])
                push!(channels, ychannel[((i, j), tick_index(i))])
            end
        end
        for i in ticking
            R[xidx[i], :] = xnew[i]
            for o in axes(ynew[i], 1)
                haskey(hidx, (i, o)) && (R[hidx[(i, o)], :] = ynew[i][o, :])
            end
        end
        if ci != 0
            Uc = continuous_inputs()
            if analysis
                ycpost = systems[ci].C * R[xidx[ci], :] + systems[ci].D * Uc
                for j in 1:partitions[ci].ny_user
                    push!(rows, ycpost[j:j, :])
                    push!(channels, ychannel[((ci, j), tick_index(ci))])
                end
            end
            Δt = (k < K ? event_times[k+1] : event_times[1] + Tc) - te
            nc = systems[ci].nx
            Φ = exp([systems[ci].A systems[ci].B; zeros(systems[ci].nu, nc + systems[ci].nu)] * Δt)
            R[xidx[ci], :] = Φ[1:nc, 1:nc] * R[xidx[ci], :] + Φ[1:nc, nc+1:end] * Uc
        end
        CD = isempty(rows) ? zeros(0, n + NU) : reduce(vcat, rows)
        push!(A, R[:, 1:n]); push!(B, R[:, n+1:end])
        push!(C, CD[:, 1:n]); push!(D, CD[:, n+1:end])
        push!(output_channels, channels)
    end
    MultirateLinearization(A, B, C, D, output_channels, Tc, t0, event_times, x_names, inports, outports, hl)
end

function multirate_linearization(sys::ModelingToolkit.AbstractSystem, inputs, outputs;
        oversample = 1, t0 = 0.0, continuous_rate = :fastest, tol = 1e-9, kwargs...)
    inputs, outputs = vcat(inputs), vcat(outputs)
    input_names = signal_channel_names(sys, inputs)
    output_names = signal_channel_names(sys, outputs)
    hl = hybrid_linearization(ModelingToolkit.linearize, sys, inputs, outputs; kwargs...)
    multirate_linearization(hl; oversample, t0, continuous_rate, tol, inputs, outputs, input_names, output_names)
end

# The matrices of the system lifted over one period (Meyer and Burrus, 1975)
function lifted_matrices(ml::MultirateLinearization)
    n = length(ml.x)
    NU = size(ml.B[1], 2)
    NY = sum(p -> length(p.channels), ml.outputs; init = 0)
    X = Matrix{Float64}(I, n, n)
    Γ = zeros(n, NU)
    H = zeros(NY, n)
    E = zeros(NY, NU)
    for k in eachindex(ml.A)
        ch = ml.output_channels[k]
        H[ch, :] = ml.C[k] * X
        E[ch, :] = ml.C[k] * Γ + ml.D[k]
        Γ = ml.A[k] * Γ + ml.B[k]
        X = ml.A[k] * X
    end
    X, Γ, H, E
end

channel_names(ports::Vector{MultiratePort}) = reduce(vcat, (length(p.times) == 1 ? [p.name] : [Symbol(p.name, "@", r) for r in eachindex(p.times)] for p in ports); init = Symbol[])

"""
    lifted_ss(ml::MultirateLinearization)

The time-invariant system that is equivalent to the periodic model `ml`, obtained by lifting
over one period: a discrete-time `NamedStateSpace` with the common period as sample
interval, in which a signal with `n` instants per period is represented by the `n` channels
`s@1, …, s@n` in chronological order. The state is the state of `ml` before the first event
of the period, and the state matrix is the monodromy matrix.

The H∞ norm of the lifted system is the induced ℓ2 norm between the sequences of the sampled
signals.
"""
function lifted_ss(ml::MultirateLinearization)
    F, G, H, E = lifted_matrices(ml)
    nsys = named_ss(ss(F, G, H, E, ml.Tc); x = ml.x, u = channel_names(ml.inputs), y = channel_names(ml.outputs))
    RobustAndOptimalControl.set_extra!(nsys, :multirate_linearization, ml)
    nsys
end

"""
    multipliers(ml::MultirateLinearization)

The characteristic multipliers of the periodic model: the eigenvalues of the monodromy
matrix. The model is exponentially stable if and only if all multipliers have magnitude
smaller than one. Hold variables that are redundant contribute multipliers that are zero.
"""
multipliers(ml::MultirateLinearization) = eigvals(lifted_matrices(ml)[1])

"""
    isstable(ml::MultirateLinearization)

Return `true` if the periodic model is exponentially stable, that is, if all characteristic
[`multipliers`](@ref) have a magnitude smaller than one.
"""
ControlSystemsBase.isstable(ml::MultirateLinearization) = all(λ -> abs(λ) < 1, multipliers(ml))

function find_port(ports, s, kind)
    s isa Integer && return ports[s]
    i = findfirst(p -> p.name == Symbol(s), ports)
    i === nothing && throw(ArgumentError("No $kind named $s, the available names are $(join((p.name for p in ports), ", "))."))
    ports[i]
end

function uniform_port(port)
    n = length(port.times)
    n == 1 && return port
    d = diff(port.times)
    all(x -> isapprox(x, d[1]; rtol = 1e-6), d) || throw(ArgumentError("The instants of the signal $(port.name) are not uniformly spaced."))
    port
end

tone(port, w) = cis.(w .* port.times)

# Frequency response of the lifted system between two signals, an array of size
# (instants of y) × (instants of u) × length(w)
function lifted_freqresp(ml, u, y, w)
    F, G, H, E = lifted_matrices(ml)
    sys = ss(F, G[:, u.channels], H[y.channels, :], E[y.channels, u.channels], ml.Tc)
    freqresp(sys, w)
end

"""
    harmonic_response(ml::MultirateLinearization, u, y, w, k)

The complex amplitude of the component with angular frequency `w + 2πk/Tc` in the
steady-state response of the output `y` to the input `u = exp(im*w*t)`, for each element of
the frequency vector `w`. `u` and `y` are names or indices of signals of `ml`. The signals
are required to have uniformly spaced instants.

`k = 0` gives the [`fundamental_response`](@ref). The components with `k ≠ 0` are caused by
sampling and holding and are absent in a time-invariant system. An output with `n` instants
per period distinguishes `n` components, `k = 0, …, n - 1`; other values of `k` coincide
with one of these.
"""
function harmonic_response(ml::MultirateLinearization, u, y, w::AbstractVector, k::Integer)
    up = uniform_port(find_port(ml.inputs, u, "input"))
    yp = uniform_port(find_port(ml.outputs, y, "output"))
    G = lifted_freqresp(ml, up, yp, w)
    Ω = 2π / ml.Tc
    map(eachslice(G; dims = 3), w) do Gi, wi
        dot(tone(yp, wi + k * Ω), Gi, tone(up, wi)) / length(yp.times)
    end
end

"""
    fundamental_response(ml::MultirateLinearization, u, y, w)

The complex amplitude of the component with the angular frequency of the input in the
steady-state response of the output `y` to the input `u = exp(im*w*t)`, for each element of
the frequency vector `w`.

For a time-invariant system this is the frequency response. For a sampled-data system it is
one entry of the harmonic response: it is not a transfer function, and neither closed-loop
stability nor stability margins can be concluded from it. The result refers to the input
held on its grid and the output observed on its grid; with a fine grid for signals of the
continuous-time partition (`oversample`), the result for a sensitivity function approaches
the fundamental sensitivity function of Goodwin and Salgado (1994).

See also [`harmonic_response`](@ref) and [`performance_frequency_gain`](@ref).
"""
fundamental_response(ml::MultirateLinearization, u, y, w::AbstractVector) = harmonic_response(ml, u, y, w, 0)

"""
    performance_frequency_gain(ml::MultirateLinearization, u, y, w)

The ratio between the root-mean-square value of the steady-state response of the output `y`
and the root-mean-square value of the input `u = exp(im*w*t)`, for each element of the
frequency vector `w`. All harmonic components of the response are included, as opposed to
[`fundamental_response`](@ref), the magnitude of which is a lower bound.
"""
function performance_frequency_gain(ml::MultirateLinearization, u, y, w::AbstractVector)
    up = uniform_port(find_port(ml.inputs, u, "input"))
    yp = uniform_port(find_port(ml.outputs, y, "output"))
    G = lifted_freqresp(ml, up, yp, w)
    map(eachslice(G; dims = 3), w) do Gi, wi
        norm(Gi * tone(up, wi)) / sqrt(length(yp.times))
    end
end

"""
    robust_frequency_gain(ml::MultirateLinearization, u, y, w)

The largest singular value of the frequency response of the lifted system between the
signals `u` and `y`, for each element of the frequency vector `w`: the largest gain, in
terms of the energy of the sequences of sampled values, over all inputs that are composed of
the frequencies `w + 2πk/Tc`. The result is periodic in `w` with period `2π/Tc`, and for a
stable model its supremum is the H∞ norm of the lifted system between the two signals, that
is, the induced ℓ2 norm between the sequences of their sampled values.
For signals with `n_u` and `n_y` instants per period,
`performance_frequency_gain ≤ sqrt(n_u/n_y)*robust_frequency_gain`.
"""
function robust_frequency_gain(ml::MultirateLinearization, u, y, w::AbstractVector)
    up = find_port(ml.inputs, u, "input")
    yp = find_port(ml.outputs, y, "output")
    G = lifted_freqresp(ml, up, yp, w)
    map(opnorm, eachslice(G; dims = 3))
end

# Search from the nominal gain towards `limit` for the first gain at which `stable` is false.
# `singular` contains gains at which the loop is not well posed; the nearest one in the
# direction of the search bounds the interval, and the search covers the gains up to it.
# Returns the boundary, or `Inf`/`0.0` if none is found before `limit`.
function stability_limit(stable, limit, singular; step = 1.002, iterations = 60)
    up = limit > 1
    candidates = filter(g -> up ? 1 < g <= limit : limit <= g < 1, singular)
    bound = isempty(candidates) ? nothing : up ? minimum(candidates) : maximum(candidates)
    # The last gain that is examined: next to the singular gain, on the side of the nominal gain
    last = bound === nothing ? limit : up ? bound * (1 - 1e-9) : bound * (1 + 1e-9)
    function bisect(lo, hi)
        for _ in 1:iterations
            mid = sqrt(lo * hi)
            stable(mid) ? (lo = mid) : (hi = mid)
        end
        lo
    end
    g = 1.0
    while true
        gn = up ? min(g * step, last) : max(g / step, last)
        stable(gn) || return bisect(g, gn)
        gn == last && return bound === nothing ? (up ? Inf : 0.0) : bound
        g = gn
    end
end

"""
    multirate_margins(S::MultirateLinearization; σ = 0, N = 400)
    multirate_margins(sys, ap; σ = 0, N = 400, kwargs...)

Stability margins of a periodic loop at a break point. `S` is the sensitivity function at
the break point, with one signal as input and output, for example created from
`get_sensitivity(sys, ap; hybrid = true)`. The second method computes the sensitivity
function at the analysis point `ap` of `sys`; `kwargs` are passed to
[`multirate_linearization`](@ref) and to the linearization.

The signal at the break point has `n` instants per period, and a perturbation of it is
classified by its action on the `n` lifted channels. The result is a named tuple with one
margin per class of perturbations:

- `gain_interval`: A numerical estimate of the interval of constant real positive gains at the break point, containing the gain one, for which the loop is stable and well posed. Gains at which the loop is not well posed are computed from the feedthrough matrix; the loss of stability is located by a search with relative step `0.002` from the gain one, followed by bisection, which can miss an interval of unstable gains that is narrower than the step. A limit of `0` or `Inf` states that no boundary was found within the search range `(1e-4, 1e4)`.
- `disk_margin_lti_lower_bound`: A `RobustAndOptimalControl.Diskmargin` for time-invariant dynamic perturbations at the rate of the signal. Such a perturbation takes values at the `n` frequencies that are aliases of each other, which are independent for a complex perturbation without further restrictions, and the margin is computed from an upper bound of the structured singular value with `n` independent complex scalars. It is a lower bound of the margin for this class on the grid of `N` frequencies; a peak between the frequencies of the grid is not detected. A perturbation with real coefficients, a fixed order, or the form of a delay belongs to a smaller class.
- `disk_margin_dynamic`: A `Diskmargin` for arbitrary dynamic perturbations of the lifted signal, with the radius `1/‖S + (σ - 1)/2‖∞`, where `S` is the lifted sensitivity function. This class includes periodically time-varying perturbations, and also perturbations that are not causal within the period. It is conservative for time-invariant perturbations.

For `n = 1` the loop at the break point is time invariant and both disk margins are computed
from the H∞ norm. Their radius agrees with `RobustAndOptimalControl.diskmargin(L, σ)` for the
loop-transfer function `L = 1/S - 1`, which is the negative of the loop-transfer function
that `ModelingToolkit.get_looptransfer` returns. `S` is required to be stable.

All margins refer to the signal as it is represented in `S`. A signal of the
continuous-time partition is represented by its values on the grid selected by
`continuous_rate` and `oversample`: the perturbation acts between the observation of the
signal on the grid and a hold on the grid. The margins then apply to the continuous-time
signal itself only if moving the perturbation to the grid changes no value that is sampled
from the signal, including the values that are read at an instant that coincides with an
instant of the grid.

`σ` is the skew of the disk.
"""
function multirate_margins(ml::MultirateLinearization; σ = 0, N = 400)
    length(ml.inputs) == 1 && length(ml.outputs) == 1 && length(ml.inputs[1].times) == length(ml.outputs[1].times) &&
        maximum(abs, ml.inputs[1].times - ml.outputs[1].times) <= 1e-9 * ml.Tc ||
        throw(ArgumentError("The sensitivity function is required to have one signal as input and output."))
    port = uniform_port(ml.inputs[1])
    n = length(port.times)
    F, G, H, E = lifted_matrices(ml)
    # With the loop-transfer function L = inv(S) - I, a real gain g at the break point gives
    # the sensitivity function S*inv(g*I + (1 - g)*S), a static feedback around S
    function stable(g)
        M = lu(g * I + (1 - g) * E; check = false)
        issuccess(M) || return false
        maximum(abs, eigvals(F - (1 - g) * G * (M \ H)); init = 0.0) < 1
    end
    stable(1.0) || throw(ArgumentError("The nominal closed loop is not stable."))
    # g*I + (1 - g)*E is singular for g = e/(e - 1), where e is a real eigenvalue of E
    singular = [real(e) / (real(e) - 1) for e in eigvals(E) if abs(imag(e)) < 1e-10 && !isapprox(real(e), 1)]
    gain_interval = (stability_limit(stable, 1e-4, singular), stability_limit(stable, 1e4, singular))

    Ω = 2π / ml.Tc
    M0 = ss(F, G, H, E + (σ - 1) / 2 * I, ml.Tc)
    w = range(0, Ω, length = N + 1)[1:end-1]
    M = freqresp(M0, w)
    # The H∞ norm, and the peak on the frequency grid as a safeguard for the norm computation
    hinf, ωinf = hinfnorm(M0)
    gridpeak, ipeak = findmax(opnorm, eachslice(M; dims = 3))
    gridpeak > hinf && ((hinf, ωinf) = (gridpeak, w[ipeak]))
    disk_margin_dynamic = Diskmargin(1 / hinf, σ; ω0 = ωinf)
    if n == 1
        disk_margin_lti_lower_bound = disk_margin_dynamic
    else
        # Change of basis from the instants within the period to the aliased frequencies,
        # W(w) = W0*D(w) with D(w) = Diagonal(cis.(-w .* times)), applied in place
        times = port.times
        W0 = [cis(-m * Ω * times[r]) / sqrt(n) for m in 0:n-1, r in 1:n]
        buffer = similar(W0)
        for (Mi, wi) in zip(eachslice(M; dims = 3), w)
            for c in 1:n, r in 1:n
                Mi[r, c] *= cis(-wi * (times[r] - times[c]))
            end
            mul!(buffer, W0, Mi)
            mul!(Mi, buffer, W0')
        end
        μ = structured_singular_value(M)
        i1 = argmax(μ)
        disk_margin_lti_lower_bound = Diskmargin(1 / μ[i1], σ; ω0 = w[i1])
    end
    (; gain_interval, disk_margin_lti_lower_bound, disk_margin_dynamic)
end

function multirate_margins(sys::ModelingToolkit.AbstractSystem, ap; σ = 0, N = 400, oversample = 1, t0 = 0.0, continuous_rate = :fastest, tol = 1e-9, kwargs...)
    hl = hybrid_linearization(get_sensitivity, sys, ap; kwargs...)
    ml = multirate_linearization(hl; oversample, t0, continuous_rate, tol)
    multirate_margins(ml; σ, N)
end

"""
    single_rate_ss(hl; Ts, method = :tustin, w_prewarp = 0)
    single_rate_ss(sys, inputs, outputs; Ts, method = :tustin, w_prewarp = 0, kwargs...)

A time-invariant approximation of a hybrid model with several clocks, obtained by converting
all partitions of the `ModelingToolkit.HybridLinearization` `hl` to one common rate and
closing the connections between them. `Ts = 0` gives a continuous-time system, in which each
clocked partition is converted with `d2c`. With `Ts > 0`, the continuous-time partition is
discretized with `c2d`, a clocked partition with the sample interval `Ts` is used as it is,
and each other clocked partition is converted to continuous time and discretized again.
`method` and `w_prewarp` are passed to `c2d` and `d2c`. The phases of the clocks are not
represented.

The result is an approximation: a model with several clocks is periodically time varying,
and aliasing between the rates is not represented. In the continuous-time result the holds
are not represented either. No stability or robustness property of the model follows from
the result. [`multirate_linearization`](@ref) provides an exact representation for
commensurate periods; `single_rate_ss` is also applicable when the periods are not
commensurate.
"""
function single_rate_ss(hl; Ts, method = :tustin, w_prewarp = 0, inputs = nothing, outputs = nothing, kwargs...)
    Ts >= 0 || throw(ArgumentError("`Ts` must be nonnegative."))
    convert_partition(P, p) = if p.Ts === nothing
        Ts == 0 ? P : c2d(P, Ts, method; w_prewarp)
    elseif Ts != 0 && isapprox(p.Ts, Ts; rtol = 1e-9)
        ss(P.A, P.B, P.C, P.D, Ts)
    else
        Pc = d2c(P, method; w_prewarp)
        Ts == 0 ? Pc : c2d(Pc, Ts, method; w_prewarp)
    end
    lsys, uoffsets, yoffsets, _ = hybrid_ss(hl; convert_partition, kwargs...)
    partitions = hl.partitions
    uinds = user_signal_indices([p.inputs[1:p.nu_user] for p in partitions], uoffsets, inputs, "inputs")
    yinds = user_signal_indices([p.outputs[1:p.ny_user] for p in partitions], yoffsets, outputs, "outputs")
    lsys[yinds, uinds]
end

function single_rate_ss(sys::ModelingToolkit.AbstractSystem, inputs, outputs; Ts, method = :tustin, w_prewarp = 0, kwargs...)
    inputs, outputs = vcat(inputs), vcat(outputs)
    hl = hybrid_linearization(ModelingToolkit.linearize, sys, inputs, outputs; kwargs...)
    lsys = single_rate_ss(hl; Ts, method, w_prewarp, inputs, outputs)
    named_ss(lsys; u = signal_channel_names(sys, inputs), y = signal_channel_names(sys, outputs), name = string(nameof(sys)))
end
