struct Path
    next_gate::Int
    pauli::PauliWord
    coeff::Float64
    priority::Float64
    seq_id::UInt64
end
"Max heap, with FIFO ordering for equal priorities."
mutable struct Frontier
    nodes::Vector{Path}
end
Frontier() = Frontier(Path[])
ahead(a::Path,b::Path) = a.priority > b.priority || (a.priority == b.priority && a.seq_id < b.seq_id)
function heappush!(h::Frontier,p::Path)
    push!(h.nodes,p)
    k = length(h.nodes)
    while k > 1
        j = k >> 1
        ahead(p,h.nodes[j]) || break
        h.nodes[k] = h.nodes[j]
        k = j
    end
    h.nodes[k] = p
end
function heappop!(h::Frontier)
    top = h.nodes[1]
    p = pop!(h.nodes)
    isempty(h.nodes) && return top
    k = 1
    while 2k <= length(h.nodes)
        j = 2k
        j < length(h.nodes) && ahead(h.nodes[j+1],h.nodes[j]) && (j += 1)
        ahead(h.nodes[j],p) || break
        h.nodes[k] = h.nodes[j]
        k = j
    end
    h.nodes[k] = p
    top
end
