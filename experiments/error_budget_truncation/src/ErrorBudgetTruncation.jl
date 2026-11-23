module ErrorBudgetTruncation

using Statistics
const B = Main.BestFirstPauli

export SelectionDiagnostics, truncate_sum, allocation_schedule, propagate_budgeted,
       fixed_threshold, expectation, policy_names

struct SelectionDiagnostics
    retained::Dict{B.PauliWord,Float64}
    discarded_mass::Float64
    discarded_terms::Int
    retained_terms::Int
    selection_seconds::Float64
    allocation::Float64
    unused_allocation::Float64
end

"Stable key used to make equal-magnitude coefficient decisions reproducible."
wordkey(p::B.PauliWord) = (p.x, p.z)

function merged_sum(terms)
    out = Dict{B.PauliWord,Float64}()
    for (p,c) in terms
        isfinite(c) || throw(ArgumentError("non-finite Pauli coefficient"))
        c == 0 && continue
        out[p] = get(out,p,0.0) + c
    end
    filter!(entry->entry.second != 0.0, out)
    out
end

"Discard a deterministic prefix of ascending |coefficient| terms within allowance."
function truncate_sum(terms::AbstractDict{B.PauliWord,<:Real}, allowance::Real)
    a = Float64(allowance)
    isfinite(a) && a >= 0 || throw(ArgumentError("allowance must be finite and nonnegative"))
    t0 = time_ns()
    merged = merged_sum(terms)
    ordered = sort(collect(merged); by=x->(abs(x[2]), wordkey(x[1])))
    retained = Dict{B.PauliWord,Float64}()
    discarded = 0.0; nd = 0
    # A term is discarded only if adding its exact mass remains inside the budget.
    # The small tolerance handles a one-ulp sum without permitting a real overrun.
    for (p,c) in ordered
        mass = abs(c)
        if discarded + mass <= a + 2eps(Float64)*max(1.0,a,discarded,mass)
            discarded += mass; nd += 1
        else
            retained[p] = c
        end
    end
    # Recompute from the selected terms for a conservative, auditable value.
    discarded = sum((abs(c) for (p,c) in ordered if !haskey(retained,p)); init=0.0)
    SelectionDiagnostics(retained, discarded, nd, length(retained),
        (time_ns()-t0)/1e9, a, max(0.0,a-discarded))
end

policy_names() = (:uniform, :front_loaded, :back_loaded, :cost_aware)

function allocation_schedule(policy::Symbol, checkpoints::Int; front_power=1.0, back_power=1.0)
    checkpoints >= 0 || throw(ArgumentError("negative checkpoint count"))
    checkpoints == 0 && return Float64[]
    policy in policy_names() || throw(ArgumentError("unknown policy $policy"))
    policy == :uniform && return fill(1.0/checkpoints, checkpoints)
    x = collect(1:checkpoints) ./ checkpoints
    w = policy == :front_loaded ? (1 .- x).^front_power :
        policy == :back_loaded ? x.^back_power : fill(1.0, checkpoints)
    # Avoid an exactly-zero first/last weight: a checkpoint must still be able to
    # spend carried-forward budget if it is the only useful checkpoint.
    w .+= eps(Float64)
    w ./ sum(w)
end

function gate_step(sum::Dict{B.PauliWord,Float64}, g)
    merged_sum((q,cq) for (p,c) in sum for (q,cq) in B.expand(g,p,c))
end

function choose_allowance(policy, remaining, k, ncheck, previous_mass, previous_terms)
    remcheck = ncheck-k+1
    if policy == :uniform
        return remaining/remcheck
    elseif policy == :front_loaded
        weights = allocation_schedule(:front_loaded, ncheck)
        return remaining * weights[k:end][1] / sum(weights[k:end])
    elseif policy == :back_loaded
        weights = allocation_schedule(:back_loaded, ncheck)
        return remaining * weights[k:end][1] / sum(weights[k:end])
    elseif policy == :cost_aware
        # Transparent online heuristic: larger current sums and recent growth get
        # more of the remaining allowance; gate count is already represented by
        # equal residual scheduling. No future circuit information is consulted.
        growth = previous_terms <= 0 ? 1.0 : max(1.0, length(previous_mass)/previous_terms)
        pressure = max(1.0, length(previous_mass)) * growth
        return remaining * pressure / (pressure + remcheck - 1)
    end
    throw(ArgumentError("unknown policy $policy"))
end

"Propagate the Heisenberg Pauli sum, merging before each certified truncation."
function propagate_budgeted(circuit, obs; epsilon_total::Real=0.0, policy::Symbol=:uniform,
        checkpoint_every::Int=1, max_terms::Int=typemax(Int), max_seconds::Real=Inf)
    eps_total = Float64(epsilon_total)
    isfinite(eps_total) && eps_total >= 0 || throw(ArgumentError("epsilon_total must be finite and nonnegative"))
    checkpoint_every >= 1 || throw(ArgumentError("checkpoint_every must be positive"))
    max_terms >= 1 || throw(ArgumentError("max_terms must be positive"))
    isfinite(max_seconds) || max_seconds == Inf || throw(ArgumentError("invalid time limit"))
    policy in policy_names() || throw(ArgumentError("unknown policy $policy"))
    B.validword(obs,circuit.nqubits)
    gates = circuit.gates
    checkpoints = cld(length(gates),checkpoint_every)
    remaining = eps_total; used = 0.0; sum = Dict(obs=>1.0); history = NamedTuple[]
    total_selection = 0.0; peak_terms = 1; status = :complete; reason = :finished
    started = time_ns(); prev_terms = 1; prev_sum = sum
    for k in eachindex(gates)
        if (time_ns()-started)/1e9 > max_seconds
            status=:incomplete; reason=:resource_limit; break
        end
        sum = gate_step(sum,gates[length(gates)-k+1])
        peak_terms = max(peak_terms,length(sum))
        if k % checkpoint_every == 0 || k == length(gates)
            checkpoint = cld(k,checkpoint_every)
            allowance = eps_total == 0 ? 0.0 : choose_allowance(policy,remaining,checkpoint,checkpoints,prev_sum,prev_terms)
            d = truncate_sum(sum,allowance)
            # The selection routine itself is the only discard point in this module.
            used += d.discarded_mass
            used <= eps_total + 16eps(Float64)*max(1.0,eps_total) || error("error ledger exceeded")
            remaining = max(0.0,eps_total-used)
            total_selection += d.selection_seconds
            sum = d.retained
            push!(history,(checkpoint=checkpoint,gate_index=k,terms_before=d.retained_terms+d.discarded_terms,
                retained_terms=d.retained_terms,discarded_terms=d.discarded_terms,
                allocated_allowance=allowance,discarded_mass=d.discarded_mass,
                error_used=used,remaining_allowance=remaining,selection_seconds=d.selection_seconds))
            prev_sum=sum; prev_terms=max(1,d.retained_terms)
            length(sum) <= max_terms || (status=:incomplete; reason=:term_limit; break)
        end
    end
    status == :complete && length(history) < checkpoints && (status=:incomplete; reason=:resource_limit)
    (;sum, expectation=expectation(sum), error_used=used, remaining_allowance=remaining,
       certified_bound=used, history, peak_terms, selection_seconds=total_selection,
       elapsed_seconds=(time_ns()-started)/1e9, status, reason,
       retained_terms=length(sum), discarded_terms=Base.sum(h.discarded_terms for h in history))
end

"Instrumented fixed-threshold baseline with the same merged representation."
function fixed_threshold(circuit, obs; min_abs_coeff::Real=0.0, max_seconds::Real=Inf)
    t=Float64(min_abs_coeff); isfinite(t) && t >= 0 || throw(ArgumentError("invalid threshold"))
    sum=Dict(obs=>1.0); discarded=0.0; hist=NamedTuple[]; started=time_ns(); status=:complete
    for k in eachindex(circuit.gates)
        sum=gate_step(sum,circuit.gates[length(circuit.gates)-k+1])
        before=length(sum); drop=0
        for (p,c) in collect(sum)
            if abs(c) < t
                discarded += abs(c); delete!(sum,p); drop += 1
            end
        end
        push!(hist,(checkpoint=k,terms_before=before,retained_terms=length(sum),discarded_terms=drop,discarded_mass=discarded))
        if (time_ns()-started)/1e9 > max_seconds; status=:incomplete; break; end
    end
    (;sum,expectation=expectation(sum),total_discarded_mass=discarded,history=hist,
      elapsed_seconds=(time_ns()-started)/1e9,status,retained_terms=length(sum))
end

expectation(terms::AbstractDict) = Base.sum((c for (p,c) in terms if p.x == 0); init=0.0)

end
