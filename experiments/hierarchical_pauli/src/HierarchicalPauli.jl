module HierarchicalPauli

using Random, LinearAlgebra, Statistics, JSON3
const B = Main.BestFirstPauli
const TermMap = Dict{B.PauliWord,Float64}
export Block, BlockMap, build_exact_block_map, compose_maps, apply_map, propagate,
       truncate_terms, strategies, recursive_map, exact_sequential_map,
       apply_map_sum, propagate_map, merge_terms,
       benchmark_record, map_equal

bitword(a,q) = !iszero(a & (UInt64(1) << (q-1)))
struct Block
    gates::Vector{B.SupportedGate}
    support::Vector{Int}
end
function Block(gates::Vector{B.SupportedGate})
    qs=Int[]
    for g in gates
        if g isa B.Rotation
            for q in 1:64
                bitword(g.generator.x | g.generator.z,q) && push!(qs,q)
            end
        else
            push!(qs,g.a); g.kind == :CNOT && push!(qs,g.b)
        end
    end
    Block(gates,sort(unique(qs)))
end
struct BlockMap
    columns::Dict{B.PauliWord,TermMap}
    support::Vector{Int}
    generated::Int
    unique_outputs::Int
    cancelled::Int
    below_threshold::Int
    compose_seconds::Float64
end
empty_map(p)=TermMap(p=>1.0)
function merge_terms(terms)
    d=TermMap()
    for (p,c) in terms
        c == 0 && continue
        d[p]=get(d,p,0.0)+c
    end
    filter!(x->x.second != 0.0,d); d
end
function gate_map(g,inputs)
    out=TermMap(); generated=0
    for (p,c) in inputs
        for (q,cq) in B.expand(g,p,c)
            generated+=1; out[q]=get(out,q,0.0)+cq
        end
    end
    merge_terms(out),generated
end
function apply_gate_sequence(gates,p)
    terms=empty_map(p); generated=0; peak=1
    for g in reverse(gates)
        terms,n=gate_map(g,terms); generated+=n; peak=max(peak,length(terms))
    end
    terms,generated,peak
end
function local_words(support)
    out=B.PauliWord[]
    for code in 0:(4^length(support)-1)
        x=z=UInt64(0); t=code
        for q in support
            v=mod(t,4); t ÷= 4; b=UInt64(1)<<(q-1)
            v==1 && (x|=b); v==2 && (z|=b); v==3 && ((x|=b);(z|=b))
        end
        push!(out,B.PauliWord(x,z))
    end
    out
end
function build_exact_block_map(gates::Vector{B.SupportedGate}; target_qubits=nothing, threshold=0.0)
    block=Block(gates); support=target_qubits === nothing ? block.support : sort(unique(target_qubits))
    cols=Dict{B.PauliWord,TermMap}(); generated=0; uniqueout=0; cancelled=0; below=0; t=time_ns()
    for p in local_words(support)
        d,n,_=apply_gate_sequence(gates,p); generated+=n; uniqueout+=length(d)
        below += count(abs(c)<threshold for c in values(d)); cols[p]=d
    end
    BlockMap(cols,support,generated,uniqueout,cancelled,below,(time_ns()-t)/1e9)
end
function compose_maps(left::BlockMap,right::BlockMap; threshold=0.0)
    t=time_ns(); cols=Dict{B.PauliWord,TermMap}(); generated=0; uniqueout=0; cancelled=0; below=0
    combined_support=sort(unique(vcat(left.support,right.support)))
    for p in local_words(combined_support)
        rcol=apply_map(right,p,1.0)
        terms=Tuple{B.PauliWord,Float64}[]
        for (mid,cr) in rcol
            lcol=apply_map(left,mid,1.0)
            for (q,c) in lcol; push!(terms,(q,cr*c)); generated+=1 end
        end
        raw=length(terms); d=merge_terms(terms); uniqueout+=length(d); cancelled+=raw-length(d)
        below += count(abs(c)<threshold for c in values(d)); cols[p]=d
    end
    BlockMap(cols,combined_support,generated,uniqueout,cancelled,below,(time_ns()-t)/1e9)
end
function apply_map(m::BlockMap,p,c=1.0)
    mask=zero(UInt64)
    for q in m.support; mask |= UInt64(1)<<(q-1) end
    localp=B.PauliWord(p.x & mask,p.z & mask)
    col=get(m.columns,localp,nothing); col === nothing && return TermMap(p=>c)
    TermMap(B.PauliWord((q.x & mask) | (p.x & ~mask),(q.z & mask) | (p.z & ~mask))=>c*v for (q,v) in col)
end
function map_equal(a::BlockMap,b::BlockMap; atol=1e-11)
    a.support==b.support || return false
    for p in union(keys(a.columns),keys(b.columns))
        da=get(a.columns,p,TermMap()); db=get(b.columns,p,TermMap()); keys(da)==keys(db) || return false
        all(isapprox(da[q],db[q];atol=atol,rtol=atol) for q in keys(da)) || return false
    end
    true
end
function truncate_terms(d::TermMap,threshold::Real;max_terms=typemax(Int))
    x=merge_terms(d); dropped=0.0; nd=0
    for (p,c) in collect(x)
        if abs(c)<threshold; dropped+=abs(c); delete!(x,p); nd+=1 end
    end
    if length(x)>max_terms
        keep=sort(collect(x);by=v->(-abs(v[2]),(v[1].x,v[1].z)))
        for (p,c) in keep[max_terms+1:end]; dropped+=abs(c); delete!(x,p); nd+=1 end
    end
    x,dropped,nd
end
function propagate(c::B.Circuit,obs;block_size=1,threshold=0.0,max_terms=typemax(Int))
    d=TermMap(obs=>1.0); peak=1; dropped=0.0; trunc_events=0; trunc_seconds=0.0; gs=c.gates
    for last in length(gs):-block_size:1
        first=max(1,last-block_size+1)
        for g in reverse(gs[first:last]); d,_=gate_map(g,d); peak=max(peak,length(d)); end
        t=time_ns(); d,dd,_=truncate_terms(d,threshold;max_terms=max_terms); trunc_seconds+=(time_ns()-t)/1e9
        dropped+=dd; trunc_events+=1
    end
    (;terms=d,expectation=sum((v for (p,v) in d if p.x==0);init=0.0),peak_terms=peak,
      dropped_mass=dropped,truncation_seconds=trunc_seconds,truncation_events=trunc_events)
end
function exact_sequential_map(gates; target_qubits=nothing)
    build_exact_block_map(gates;target_qubits=target_qubits)
end
function recursive_map(gates; split=:balanced, max_support=typemax(Int), max_columns=typemax(Int))
    length(gates)>0 || return build_exact_block_map(B.SupportedGate[];target_qubits=Int[])
    if length(gates)==1
        return build_exact_block_map(gates)
    end
    k=split==:left_deep ? length(gates)-1 : length(gates)÷2
    left=recursive_map(gates[1:k];split=split,max_support=max_support,max_columns=max_columns)
    right=recursive_map(gates[k+1:end];split=split,max_support=max_support,max_columns=max_columns)
    support=sort(unique(vcat(left.support,right.support)))
    if length(support)>max_support || 4^length(support)>max_columns
        # Safety fallback: retain separate blocks; caller can apply them sequentially.
        return build_exact_block_map(gates;target_qubits=support)
    end
    compose_maps(left,right)
end
function strategy_map(gates;strategy=:balanced,max_support=typemax(Int),max_columns=typemax(Int))
    strategy==:sequential && return exact_sequential_map(gates)
    recursive_map(gates;split=strategy,max_support=max_support,max_columns=max_columns)
end
function apply_map_sum(m::BlockMap,d::TermMap)
    out=TermMap()
    for (p,c) in d
        for (q,v) in apply_map(m,p,c); out[q]=get(out,q,0.0)+v end
    end
    merge_terms(out)
end
function propagate_map(c::B.Circuit,obs;strategy=:balanced,block_size=length(c.gates),threshold=0.0,
        max_support=typemax(Int),max_columns=typemax(Int))
    d=TermMap(obs=>1.0); peak=1; maps=0; compose_seconds=0.0; trunc_events=0; dropped=0.0
    gs=c.gates
    for last in length(gs):-block_size:1
        first=max(1,last-block_size+1); m=strategy_map(gs[first:last];strategy=strategy,max_support=max_support,max_columns=max_columns)
        maps+=1; compose_seconds+=m.compose_seconds; d=apply_map_sum(m,d); peak=max(peak,length(d))
        if threshold>0
            d,dd,_=truncate_terms(d,threshold); dropped+=dd; trunc_events+=1
        end
    end
    (;terms=d,expectation=sum((v for (p,v) in d if p.x==0);init=0.0),peak_terms=peak,map_count=maps,
      composition_seconds=compose_seconds,truncation_events=trunc_events,dropped_mass=dropped)
end
strategies()=(:left_deep,:balanced,:support_constrained)
make_circuit(n,family,layers,seed)=B.make_circuit(n,family,layers,seed)
observable(n,k=:Z_center)=B.observable(n,k)
dense_expectation(c,p)=B.dense_expectation(c,p)
benchmark_record(;kwargs...)=(;kwargs...,julia_version=string(VERSION),threads=Threads.nthreads(),package_version="PauliPropagation 0.9.2")
end
