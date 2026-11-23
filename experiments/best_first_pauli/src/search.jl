# Compensated accumulation for terminal contributions and frontier mass.
mutable struct Accumulator
    sum::Float64
    correction::Float64
    movement::Float64
end
Accumulator() = Accumulator(0.0,0.0,0.0)
function add!(a::Accumulator,x::Float64)
    t = a.sum + x
    a.correction += abs(a.sum) >= abs(x) ? (a.sum-t)+x : (x-t)+a.sum
    a.sum = t
    a.movement += abs(x)
end
value(a::Accumulator) = a.sum + a.correction
priority_function(mode::Symbol) = mode == :abs_coeff ? abs : throw(ArgumentError("unsupported priority_mode $mode"))
priority_function(f::Function) = f
function checked_priority(f,c)
    p = Float64(f(c))
    isfinite(p) && p >= 0 || throw(ArgumentError("priority must be finite and nonnegative"))
    p
end
"""
    bestfirst(circuit, observable; max_expansions=2_000_000, max_frontier_nodes=250_000,
              max_runtime_seconds=30.0, error_tolerance=nothing, priority_mode=:abs_coeff)

Fixed-angle |0> expectation. All gates are in Schrodinger order. A memory stop
leaves the parent intact; no child/subtree is discarded. The frontier bound is
independent of the ordering heuristic. Caps include the root; expansions exclude
terminal pops. Time is checked between operations, so one operation may overrun.
"""
function bestfirst(circuit::Circuit,obs::PauliWord; coefficient::Real=1.0,
                   max_expansions::Int=2_000_000,max_frontier_nodes::Int=250_000,
                   max_runtime_seconds::Real=30.0,error_tolerance=nothing,
                   priority_mode=:abs_coeff,measure_bookkeeping::Bool=false)
    validword(obs,circuit.nqubits)
    max_expansions >= 0 && max_frontier_nodes >= 1 || throw(ArgumentError("invalid node limits"))
    (isfinite(max_runtime_seconds) || max_runtime_seconds == Inf) && max_runtime_seconds >= 0 || throw(ArgumentError("invalid time limit"))
    error_tolerance === nothing || (isfinite(error_tolerance) && error_tolerance >= 0) || throw(ArgumentError("invalid tolerance"))
    c = Float64(coefficient)
    isfinite(c) || throw(ArgumentError("coefficient must be finite"))
    priority = priority_function(priority_mode)
    h = Frontier()
    mass,total = Accumulator(),Accumulator()
    seq = UInt64(0)
    if c != 0
        heappush!(h,Path(length(circuit.gates),obs,c,checked_priority(priority,c),seq))
        add!(mass,abs(c))
    end
    expanded = generated = leaves = pops = depth = 0
    pushes = peak = length(h.nodes)
    heap_total = heap_samples = 0
    bookkeeping_ns = kernel_ns = UInt64(0)
    status,reason = :exact,:frontier_empty
    start = time_ns()
    bound() = max(0.0,value(mass)) + 8eps(Float64)*mass.movement
    while !isempty(h.nodes)
        if error_tolerance !== nothing && bound() <= error_tolerance
            status,reason = :certified,:error_tolerance
            break
        end
        if (time_ns()-start)/1e9 >= max_runtime_seconds
            status,reason = :resource_limited,:runtime
            break
        end
        p = h.nodes[1]
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
            if length(h.nodes)-1+countchildren > max_frontier_nodes
                status,reason = :resource_limited,:frontier
                break
            end
        end
        if measure_bookkeeping
            t=time_ns(); heappop!(h); pops += 1; add!(mass,-abs(p.coeff)); bookkeeping_ns += time_ns()-t
        else
            heappop!(h); pops += 1; add!(mass,-abs(p.coeff))
        end
        depth = max(depth,length(circuit.gates)-p.next_gate)
        if p.next_gate == 0
            leaves += 1
            iszero(p.pauli.x) && add!(total,p.coeff)
        else
            expanded += 1
            t=measure_bookkeeping ? time_ns() : UInt64(0)
            for (q,cq) in children
                cq == 0 && continue
                seq += 1
                heappush!(h,Path(p.next_gate-1,q,cq,checked_priority(priority,cq),seq))
                add!(mass,abs(cq)); generated += 1; pushes += 1
            end
            measure_bookkeeping && (bookkeeping_ns += time_ns()-t)
        end
        peak = max(peak,length(h.nodes))
        heap_total += length(h.nodes); heap_samples += 1
    end
    (; expectation=value(total),status,reason,frontier_bound=isempty(h.nodes) ? 0.0 : bound(),
       frontier_abs_coeff_sum=isempty(h.nodes) ? 0.0 : max(0.0,value(mass)),
       elapsed_seconds=(time_ns()-start)/1e9,expanded_nodes=expanded,generated_children=generated,
       completed_leaves=leaves,queue_pushes=pushes,queue_pops=pops,peak_frontier=peak,
       frontier_nodes=length(h.nodes),approx_node_bytes=peak*sizeof(Path),max_tree_depth=depth,
       mean_frontier=heap_samples == 0 ? Float64(peak) : heap_total/heap_samples,
       stack_pushes=missing,stack_pops=missing,
       bookkeeping_seconds=measure_bookkeeping ? bookkeeping_ns/1e9 : missing,
       kernel_seconds=measure_bookkeeping ? kernel_ns/1e9 : missing)
end
