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
states of the clocked partitions, one hold variable per clocked signal that is read by
another partition, and one hold variable per input of the continuous-time partition.

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
    T1 * lcm(numerator.(ratios)...) / gcd(denominator.(ratios)...)
end

# The instants `phase + m*T` in `[t0, t0 + Tc)`
function tick_times(T, phase, t0, Tc; tol)
    n = round(Int, Tc / T)
    m0 = ceil(Int, (t0 - phase) / T - tol)
    [phase + (m0 + m) * T for m in 0:n-1]
end

# Order in which the clocked partitions `ticking` are evaluated at an instant: a partition
# that reads the output of another partition ticking at the same instant is evaluated after it
function evaluation_order(ticking, connections, ci)
    remaining = collect(ticking)
    order = Int[]
    while !isempty(remaining)
        i = findfirst(remaining) do p
            !any(c -> c.to == p && c.from != ci && c.from in remaining && c.from != p, connections)
        end
        i === nothing && throw(ArgumentError("The clocked partitions $(remaining) that tick at the same instant depend on each other in a cycle, the model is not causal."))
        push!(order, remaining[i])
        deleteat!(remaining, i)
    end
    order
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

The semantics at an instant where clocks tick are:
- A sampled signal of the continuous-time partition is read before any hold is updated at
  that instant.
- A clocked partition that reads a signal of another clocked partition ticking at the same
  instant reads the new value.
- All holds are updated, after which the continuous-time partition is integrated exactly to
  the next instant with constant inputs.

# Keyword arguments
- `continuous_rate`: Inputs and outputs of the continuous-time partition are held and observed on a uniform grid. With `:fastest`, the interval of the grid is the smallest clock period divided by `oversample`; with `:period`, it is the common period.
- `oversample`: An integer that refines the grid of the signals of the continuous-time partition, which makes the behavior between the ticks of the clocks visible. The ticks of the clocks are not changed.
- `t0`: The instant at which the period starts. The characteristic multipliers do not depend on `t0`, the input-output behavior does.
- `inputs`, `outputs`: The variables that were used as inputs and outputs of the linearization. They are required to recover the order of the signals if these belong to several clock partitions.
- `tol`: Relative tolerance for the recognition of commensurate periods and of coinciding instants.

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
    ticks = Dict(i => tick_times(partitions[i].Ts, partitions[i].clock.phase, t0, Tc; tol) for i in clocked)
    if ci != 0 && (partitions[ci].nu_user > 0 || partitions[ci].ny_user > 0)
        Ta = continuous_rate === :period ? Tc : minimum(periods) / oversample
        ticks[ci] = tick_times(Ta, 0.0, t0, Tc; tol)
    elseif ci != 0
        ticks[ci] = Float64[]
    end
    event_times = Float64[]
    for t in sort!(reduce(vcat, values(ticks)))
        (isempty(event_times) || t - event_times[end] > atol) && push!(event_times, t)
    end
    K = length(event_times)
    event_of(t) = findfirst(te -> abs(te - t) <= atol, event_times)

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
        ticking = [i for i in clocked if any(t -> abs(t - te) <= atol, ticks[i])]
        tick_index(i) = findfirst(t -> abs(t - te) <= atol, ticks[i])
        analysis = ci != 0 && tick_index(ci) !== nothing
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
        ynew = Dict{Int, Matrix{Float64}}()
        xnew = Dict{Int, Matrix{Float64}}()
        rows = Matrix{Float64}[]
        channels = Int[]
        for i in evaluation_order(ticking, hl.connections, ci)
            s = systems[i]
            p = partitions[i]
            Ui = zeros(s.nu, n + NU)
            for j in 1:s.nu
                if j <= p.nu_user
                    Ui[j, :] = unit(uchannel[((i, j), tick_index(i))])
                elseif haskey(source, (i, j))
                    from, o = source[(i, j)]
                    Ui[j, :] = from == ci ? yc[o, :] : haskey(ynew, from) ? ynew[from][o, :] : R[hidx[(from, o)], :]
                end
            end
            ynew[i] = s.C * R[xidx[i], :] + s.D * Ui
            xnew[i] = s.A * R[xidx[i], :] + s.B * Ui
            for j in 1:p.ny_user
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

The H∞ norm of the lifted system is the induced norm between the sequences of the sampled
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
smaller than one. Hold variables contribute multipliers that are zero.
"""
multipliers(ml::MultirateLinearization) = eigvals(lifted_matrices(ml)[1])

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
sampling and holding and are absent in a time-invariant system.
"""
function harmonic_response(ml::MultirateLinearization, u, y, w::AbstractVector, k::Integer)
    up = uniform_port(find_port(ml.inputs, u, "input"))
    yp = uniform_port(find_port(ml.outputs, y, "output"))
    G = lifted_freqresp(ml, up, yp, w)
    Ω = 2π / ml.Tc
    [dot(tone(yp, w[i] + k * Ω), G[:, :, i] * tone(up, w[i])) / length(yp.times) for i in eachindex(w)]
end

"""
    fundamental_response(ml::MultirateLinearization, u, y, w)

The complex amplitude of the component with the angular frequency of the input in the
steady-state response of the output `y` to the input `u = exp(im*w*t)`, for each element of
the frequency vector `w`.

For a time-invariant system this is the frequency response. For a sampled-data system it is
one entry of the harmonic response: it is not a transfer function, and neither closed-loop
stability nor stability margins can be concluded from it. If `ml` is a sensitivity function,
the result is the fundamental sensitivity function of Goodwin and Salgado (1994).

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
    [norm(G[:, :, i] * tone(up, w[i])) / sqrt(length(yp.times)) for i in eachindex(w)]
end

"""
    robust_frequency_gain(ml::MultirateLinearization, u, y, w)

The largest singular value of the frequency response of the lifted system between the
signals `u` and `y`, for each element of the frequency vector `w`: the largest gain over all
inputs that are composed of the frequencies `w + 2πk/Tc`. The result is periodic in `w` with
period `2π/Tc`, and its supremum is the induced norm between the sequences of the sampled
signals.
"""
function robust_frequency_gain(ml::MultirateLinearization, u, y, w::AbstractVector)
    up = find_port(ml.inputs, u, "input")
    yp = find_port(ml.outputs, y, "output")
    G = lifted_freqresp(ml, up, yp, w)
    [opnorm(G[:, :, i]) for i in eachindex(w)]
end

# Search from the nominal gain towards `limit` for the first gain at which `stable` is false
function stability_limit(stable, limit; step = 1.01, iterations = 60)
    up = limit > 1
    g = 1.0
    while true
        gn = up ? g * step : g / step
        if up ? gn > limit : gn < limit
            return up ? Inf : 0.0
        end
        if !stable(gn)
            lo, hi = g, gn
            for _ in 1:iterations
                mid = sqrt(lo * hi)
                stable(mid) ? (lo = mid) : (hi = mid)
            end
            return lo
        end
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

- `gain_interval`: The interval of constant real gains at the break point, containing the gain one, for which the loop is stable. Exact, within the search range `(1e-4, 1e4)`.
- `disk_margin_lti_lower_bound`: A `RobustAndOptimalControl.Diskmargin` for time-invariant dynamic perturbations at the rate of the signal. Such a perturbation takes independent values at the `n` frequencies that are aliases of each other, and the margin is computed from an upper bound of the structured singular value with `n` independent complex scalars. It is a lower bound of the margin on the frequency grid.
- `disk_margin_dynamic`: A `Diskmargin` for arbitrary dynamic perturbations of the lifted signal, which include periodically time-varying perturbations. It is conservative for time-invariant perturbations.

For `n = 1` the loop at the break point is time invariant and the two disk margins coincide
with `RobustAndOptimalControl.diskmargin`.

A signal of the continuous-time partition is represented by its values on the grid selected
by `continuous_rate` and `oversample`, and the perturbations act on this representation.
Only the gain interval then refers to a perturbation of the continuous-time signal itself.

`σ` is the skew of the disk and `N` the number of frequencies in the grid.
"""
function multirate_margins(ml::MultirateLinearization; σ = 0, N = 400)
    length(ml.inputs) == 1 && length(ml.outputs) == 1 && ml.inputs[1].times ≈ ml.outputs[1].times ||
        throw(ArgumentError("The sensitivity function is required to have one signal as input and output."))
    port = uniform_port(ml.inputs[1])
    n = length(port.times)
    F, G, H, E = lifted_matrices(ml)
    # With the loop-transfer function L = inv(S) - I, a real gain g at the break point gives
    # the sensitivity function S*inv(g*I + (1 - g)*S), a static feedback around S
    function stable(g)
        Q = g * I + (1 - g) * E
        M = lu(Q; check = false)
        (issuccess(M) && cond(Q) < 1e12) || return false
        maximum(abs, eigvals(F - (1 - g) * G * (M \ H))) < 1
    end
    stable(1.0) || throw(ArgumentError("The nominal closed loop is not stable."))
    gain_interval = (stability_limit(stable, 1e-4), stability_limit(stable, 1e4))

    Ω = 2π / ml.Tc
    w = range(0, Ω, length = N + 1)[1:end-1]
    S = freqresp(ss(F, G, H, E, ml.Tc), w)
    M = similar(S)
    for i in eachindex(w)
        # Change of basis from the instants within the period to the aliased frequencies
        W = [cis(-(w[i] + m * Ω) * port.times[r]) / sqrt(n) for m in 0:n-1, r in 1:n]
        M[:, :, i] = W * S[:, :, i] * W' + (σ - 1) / 2 * I
    end
    μ = n == 1 ? vec(abs.(M)) : structured_singular_value(M)
    full = [opnorm(M[:, :, i]) for i in eachindex(w)]
    i1, i2 = argmax(μ), argmax(full)
    (; gain_interval,
        disk_margin_lti_lower_bound = Diskmargin(1 / μ[i1], σ; ω0 = w[i1]),
        disk_margin_dynamic = Diskmargin(1 / full[i2], σ; ω0 = w[i2]))
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
discretized with `c2d` and each clocked partition with a different sample interval is
converted to continuous time and discretized again. `method` and `w_prewarp` are passed to
`c2d` and `d2c`.

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
