"""
    depthfirst(circuit, observable; max_expansions=2_000_000,
               max_frontier_nodes=250_000, max_runtime_seconds=30.0,
               error_tolerance=nothing)

Unmerged depth-first traversal using the same `Path`, propagation kernel,
terminal overlap, compensated sums, limits, and frontier error bound as
`bestfirst`. Children returned by `expand` are ordered cosine/unchanged first and
sine branch second. DFS pushes them in reverse order so the cosine/unchanged child
is visited first. This convention is deterministic.

Before an expansion that would exceed the frontier cap, the parent remains on the
stack. Thus every unexplored subtree remains represented in the certified bound.
"""
function depthfirst(circuit::Circuit,obs::PauliWord; coefficient::Real=1.0,
                    max_expansions::Int=2_000_000,max_frontier_nodes::Int=250_000,
                    max_runtime_seconds::Real=30.0,error_tolerance=nothing,
                    measure_bookkeeping::Bool=false)
    validword(obs,circuit.nqubits)
    max_expansions >= 0 && max_frontier_nodes >= 1 || throw(ArgumentError("invalid node limits"))
    isfinite(max_runtime_seconds) || max_runtime_seconds == Inf || throw(ArgumentError("invalid time limit"))
    max_runtime_seconds >= 0 || throw(ArgumentError("invalid time limit"))
    error_tolerance === nothing || (isfinite(error_tolerance) && error_tolerance >= 0) || throw(ArgumentError("invalid tolerance"))
    c = Float64(coefficient)
    isfinite(c) || throw(ArgumentError("coefficient must be finite"))
    stack = Path[]
    mass,total = Accumulator(),Accumulator()
    if c != 0
        push!(stack,Path(length(circuit.gates),obs,c,abs(c),0))
        add!(mass,abs(c))
    end
    expanded = generated = leaves = pops = depth = 0
    pushes = peak = length(stack)
    frontier_total = frontier_samples = 0
    bookkeeping_ns = kernel_ns = UInt64(0)
    status,reason = :exact,:frontier_empty
    start = time_ns()
    bound() = max(0.0,value(mass)) + 8eps(Float64)*mass.movement
    while !isempty(stack)
        if error_tolerance !== nothing && bound() <= error_tolerance
            status,reason = :certified,:error_tolerance
            break
        end
        if (time_ns()-start)/1e9 >= max_runtime_seconds
            status,reason = :resource_limited,:runtime
            break
        end
        p = stack[end]
        children = nothing
        if p.next_gate > 0
            if expanded >= max_expansions
                status,reason = :resource_limited,:expansions
                break
            end
            if measure_bookkeeping
                t=time_ns(); children=expand(circuit.gates[p.next_gate],p.pauli,p.coeff); kernel_ns += time_ns()-t
            else
                children=expand(circuit.gates[p.next_gate],p.pauli,p.coeff)
            end
            countchildren = count(t -> t[2] != 0,children)
            if length(stack)-1+countchildren > max_frontier_nodes
                status,reason = :resource_limited,:frontier
                break
            end
        end
        if measure_bookkeeping
            t=time_ns(); pop!(stack); pops += 1; add!(mass,-abs(p.coeff)); bookkeeping_ns += time_ns()-t
        else
            pop!(stack); pops += 1; add!(mass,-abs(p.coeff))
        end
        depth = max(depth,length(circuit.gates)-p.next_gate)
        if p.next_gate == 0
            leaves += 1
            iszero(p.pauli.x) && add!(total,p.coeff)
        else
            expanded += 1
            if measure_bookkeeping
                t=time_ns()
                for (q,cq) in Iterators.reverse(children)
                    cq == 0 && continue
                    push!(stack,Path(p.next_gate-1,q,cq,abs(cq),UInt64(pushes)))
                    add!(mass,abs(cq)); generated += 1; pushes += 1
                end
                bookkeeping_ns += time_ns()-t
            else
                for (q,cq) in Iterators.reverse(children)
                    cq == 0 && continue
                    push!(stack,Path(p.next_gate-1,q,cq,abs(cq),UInt64(pushes)))
                    add!(mass,abs(cq)); generated += 1; pushes += 1
                end
            end
        end
        peak = max(peak,length(stack))
        frontier_total += length(stack); frontier_samples += 1
    end
    elapsed = (time_ns()-start)/1e9
    (; expectation=value(total),status,reason,frontier_bound=isempty(stack) ? 0.0 : bound(),
       frontier_abs_coeff_sum=isempty(stack) ? 0.0 : max(0.0,value(mass)),
       elapsed_seconds=elapsed,expanded_nodes=expanded,generated_children=generated,
       completed_leaves=leaves,stack_pushes=pushes,stack_pops=pops,queue_pushes=missing,
       queue_pops=missing,peak_frontier=peak,frontier_nodes=length(stack),
       approx_node_bytes=peak*sizeof(Path),max_tree_depth=depth,
       mean_frontier=frontier_samples == 0 ? Float64(peak) : frontier_total/frontier_samples,
       bookkeeping_seconds=measure_bookkeeping ? bookkeeping_ns/1e9 : missing,
       kernel_seconds=measure_bookkeeping ? kernel_ns/1e9 : missing)
end
