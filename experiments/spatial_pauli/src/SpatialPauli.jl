module SpatialPauli

const B = let path = normpath(joinpath(@__DIR__, "..", "..", "best_first_pauli", "src", "BestFirstPauli.jl"))
    if !isdefined(Main, :BestFirstPauli)
        Base.include(Main, path)
    end
    Base.invokelatest(getproperty,Main,:BestFirstPauli)
end
const Word = B.PauliWord
const Factor = Dict{Word,Float64}

export spatial, sequential, expand_global, expectation_zero, partition, ResourceLimit,
       B, Component, Operator, Limits, stats, factorized, local_add, local_scale

struct ResourceLimit <: Exception
    reason::String
end
Base.showerror(io::IO, e::ResourceLimit) = print(io, e.reason)
Base.@kwdef struct Limits
    max_components::Int = 100_000
    max_local_terms::Int = 1_000_000
    max_contributions::Int = 1_000_000
    max_aggregation_passes::Int = 4
    max_representation_bytes::Int = typemax(Int)
end
struct Node
    lo::Int
    hi::Int
    left::Union{Nothing,Node}
    right::Union{Nothing,Node}
end
function tree(nodes::Vector{Node},lo::Int,hi::Int)
    lo==hi && return nodes[lo]
    mid=(lo+hi) ÷ 2
    left=tree(nodes,lo,mid); right=tree(nodes,mid+1,hi)
    Node(left.lo,right.hi,left,right)
end
function leaves(node::Node)
    node.left === nothing && return [node]
    vcat(leaves(node.left),leaves(node.right))
end
function lowest_common_ancestor(node::Node, support::UInt64)
    node.left === nothing && return node
    if (support & ~mask64(node.left))==0
        return lowest_common_ancestor(node.left,support)
    elseif (support & ~mask64(node.right))==0
        return lowest_common_ancestor(node.right,support)
    end
    node
end
function partition(n::Int,leaf_size::Int)
    1 <= n <= 64 && leaf_size >= 1 || throw(ArgumentError("invalid partition"))
    nodes=[Node(q,min(q+leaf_size-1,n),nothing,nothing) for q in 1:leaf_size:n]
    tree(nodes,1,length(nodes))
end
mask(node::Node) = ((UInt64(1) << (node.hi-node.lo+1))-1) << (node.lo-1) # special case 64 below
mask64(node::Node) = node.hi-node.lo+1 == 64 ? typemax(UInt64) : mask(node)
struct Component
    coeff::Float64
    factors::Vector{Factor}
end
function local_add(a::Factor,b::Factor)
    out=copy(a)
    for (p,v) in b; out[p]=get(out,p,0.0)+v; end
    Factor(p=>v for (p,v) in out if v!=0.0)
end
function local_scale(a::Real,f::Factor)
    out=Factor()
    for (p,v) in f
        z=Float64(a)*v
        a!=0 && v!=0 && z==0 && throw(ArgumentError("local scaling underflow"))
        z!=0 && (out[p]=z)
    end
    out
end
mutable struct Operator
    n::Int
    root::Node
    blocks::Vector{Node}
    components::Vector{Component}
    cross_gates::Int
    peak_components::Int
    peak_terms::Int
    peak_contributions::Int
    before_aggregation::Int
    aggregation_ns::Int
    limits::Limits
end
function Operator(n::Int, leaf_size::Int, terms; limits=Limits())
    1 <= n <= 64 || throw(ArgumentError("n must be in 1:64"))
    root=partition(n,leaf_size); blocks=leaves(root)
    cs=Component[]
    for (p,c) in terms
        B.validword(p,n)
        push!(cs,Component(Float64(c),[Factor(Word(p.x & mask64(b),p.z & mask64(b))=>1.0) for b in blocks]))
    end
    op=Operator(n,root,blocks,cs,0,0,0,0,0,0,limits)
    aggregate!(op)
end
function factorized(n::Int,leaf_size::Int,components::Vector{Component};limits=Limits())
    1 <= n <= 64 || throw(ArgumentError("n must be in 1:64"))
    root=partition(n,leaf_size); blocks=leaves(root)
    for c in components
        length(c.factors)==length(blocks) || throw(ArgumentError("factor count does not match leaves"))
        for (f,b) in zip(c.factors,blocks), p in keys(f)
            ((p.x | p.z) & ~mask64(b))==0 || throw(ArgumentError("factor exceeds leaf support"))
        end
    end
    op=Operator(n,root,blocks,copy(components),0,0,0,0,0,0,limits)
    aggregate!(op)
end
key(f::Factor) = Tuple(sort!(collect(f); by=x->(x.first.x,x.first.z)))
function normalize(c::Component)
    isfinite(c.coeff) || throw(ArgumentError("nonfinite component coefficient"))
    scale=c.coeff; fs=Factor[]
    for f in c.factors
        all(isfinite,values(f)) || throw(ArgumentError("nonfinite local coefficient"))
        clean=Factor(p=>v for (p,v) in f if v != 0.0)
        isempty(clean) && return nothing
        pivotword=first(keys(clean)); pivotabs=abs(clean[pivotword])
        for p in keys(clean)
            a=abs(clean[p])
            if a>pivotabs || (a==pivotabs && (p.x,p.z)<(pivotword.x,pivotword.z))
                pivotword=p; pivotabs=a
            end
        end
        pivot=clean[pivotword]
        newscale=scale*pivot
        scale!=0 && pivot!=0 && newscale==0 && throw(ArgumentError("normalization underflow"))
        scale=newscale
        isfinite(scale) || throw(ArgumentError("normalization overflow"))
        normalized=Factor()
        for (p,v) in clean
            z=v/pivot
            v!=0 && z==0 && throw(ArgumentError("normalization quotient underflow"))
            normalized[p]=z
        end
        push!(fs,normalized)
    end
    scale == 0.0 ? nothing : Component(scale,fs)
end
function aggregate!(op::Operator)
    t=time_ns(); op.before_aggregation += length(op.components)
    merged=Dict{Tuple,Component}()
    for raw in op.components
        c=normalize(raw); c === nothing && continue
        k=Tuple(key(f) for f in c.factors)
        if haskey(merged,k)
            old=merged[k]; merged[k]=Component(old.coeff+c.coeff,old.factors)
        else
            merged[k]=c
        end
    end
    op.components=[c for c in values(merged) if c.coeff != 0.0]
    # Merge along every factor position. Each pass is exact but rank minimization is not attempted.
    for _ in 1:op.limits.max_aggregation_passes
        changed=false
        for position in eachindex(op.blocks)
            groups=Dict{Tuple,Component}()
            for c in op.components
                k=Tuple(key(c.factors[i]) for i in eachindex(c.factors) if i!=position)
                if haskey(groups,k)
                    old=groups[k]; fs=copy(old.factors)
                    fs[position]=local_add(fs[position],local_scale(c.coeff,c.factors[position]))
                    groups[k]=Component(1.0,fs); changed=true
                else
                    fs=copy(c.factors); fs[position]=local_scale(c.coeff,fs[position])
                    groups[k]=Component(1.0,fs)
                end
            end
            if changed
                op.components=Component[]
                for c in values(groups)
                    v=normalize(c); v === nothing || push!(op.components,v)
                end
            end
        end
        changed || break
    end
    op.peak_components=max(op.peak_components,length(op.components))
    nt=sum(sum(length(f) for f in c.factors) for c in op.components; init=0)
    op.peak_terms=max(op.peak_terms,nt)
    length(op.components)>op.limits.max_components && throw(ResourceLimit("factorized component limit exceeded"))
    nt>op.limits.max_local_terms && throw(ResourceLimit("local term limit exceeded"))
    op.limits.max_representation_bytes!=typemax(Int) && Base.summarysize(op.components)>op.limits.max_representation_bytes && throw(ResourceLimit("representation byte estimate exceeded"))
    op.aggregation_ns += time_ns()-t
    op
end
function finalize_components!(op::Operator)
    if length(op.components)==1
        c=only(op.components)
        isfinite(c.coeff) || throw(ArgumentError("nonfinite component coefficient"))
        nt=0
        for f in c.factors
            for p in collect(keys(f))
                v=f[p]
                isfinite(v) || throw(ArgumentError("nonfinite local coefficient"))
                v==0.0 && delete!(f,p)
            end
            if isempty(f) || c.coeff==0.0
                empty!(op.components); break
            end
            nt+=length(f)
        end
        op.peak_components=max(op.peak_components,length(op.components))
        op.peak_terms=max(op.peak_terms,nt)
        nt>op.limits.max_local_terms && throw(ResourceLimit("local term limit exceeded"))
        op.limits.max_representation_bytes!=typemax(Int) && Base.summarysize(op.components)>op.limits.max_representation_bytes && throw(ResourceLimit("representation byte estimate exceeded"))
        return op
    end
    aggregate!(op)
end
function append_checked!(op::Operator, out::Vector{Component}, c::Component)
    push!(out,c)
    length(out)>op.limits.max_contributions && throw(ResourceLimit("intermediate contribution limit exceeded"))
end
function support(g::B.Rotation)
    g.generator.x | g.generator.z
end
support(g::B.Clifford) = (UInt64(1) << (g.a-1)) | (g.kind==:CNOT ? UInt64(1) << (g.b-1) : UInt64(0))
function affected(op::Operator,g)
    s=support(g)
    findall(b->(s & mask64(b)) != 0,op.blocks)
end
function local_gate(op::Operator,g,i)
    out=Component[]
    for c in op.components
        f=Factor()
        for (p,v) in c.factors[i], (q,w) in B.expand(g,p,v)
            w==0.0 && continue
            f[q]=get(f,q,0.0)+w
        end
        fs=copy(c.factors); fs[i]=f
        append_checked!(op,out,Component(c.coeff,fs))
    end
    op.components=out
end
function split_parity(f::Factor, gen::Word)
    a=(Factor(),Factor())
    for (p,v) in f
        i=B.commutes(gen,p) ? 1 : 2
        a[i][p]=v
    end
    a
end
function rotated_factor(f::Factor, gen::Word, anti::Bool)
    out=Factor()
    for (p,v) in f
        q,phase=B.product(gen,p)
        s=anti ? imag(phase) : real(phase)
        out[q]=get(out,q,0.0)+v*s
    end
    out
end
function cross_rotation!(op::Operator,g::B.Rotation,ids)
    out=Component[]
    for c in op.components
        branches=[(copy(c.factors),0)]
        for i in ids
            gen=Word(g.generator.x & mask64(op.blocks[i]),g.generator.z & mask64(op.blocks[i]))
            next=Tuple{Vector{Factor},Int}[]
            for (fs,k) in branches
                pair=split_parity(fs[i],gen)
                for anti in 0:1
                    isempty(pair[anti+1]) && continue
                    copyfs=copy(fs); copyfs[i]=pair[anti+1]
                    push!(next,(copyfs,k+anti))
                    length(next)>op.limits.max_contributions && throw(ResourceLimit("parity split limit exceeded"))
                end
            end
            branches=next
        end
        for (fs,k) in branches
            if iseven(k)
                append_checked!(op,out,Component(c.coeff,fs))
            else
                g.cosine != 0 && append_checked!(op,out,Component(c.coeff*g.cosine,fs))
                transformed=copy(fs)
                for i in ids
                    gen=Word(g.generator.x & mask64(op.blocks[i]),g.generator.z & mask64(op.blocks[i]))
                    anti=!B.commutes(gen,first(keys(fs[i])))
                    transformed[i]=rotated_factor(fs[i],gen,anti)
                end
                sign=isodd((k+1)÷2) ? -1.0 : 1.0
                g.sine != 0 && append_checked!(op,out,Component(c.coeff*g.sine*sign,transformed))
            end
        end
    end
    op.components=out
end
function letter(p::Word,q::Int)
    x=B.bit(p.x,q); z=B.bit(p.z,q)
    x ? (z ? :Y : :X) : (z ? :Z : :I)
end
function setletter(p::Word,q::Int,s::Symbol)
    bit=UInt64(1) << (q-1)
    Word((p.x & ~bit) | (s in (:X,:Y) ? bit : 0), (p.z & ~bit) | (s in (:Z,:Y) ? bit : 0))
end
function cross_cnot!(op::Operator,g::B.Clifford,ids)
    ia=findfirst(b->b.lo<=g.a<=b.hi,op.blocks)
    ib=findfirst(b->b.lo<=g.b<=b.hi,op.blocks)
    (ia===nothing || ib===nothing) && throw(ArgumentError("CNOT outside partition"))
    ia==ib && throw(ArgumentError("cross CNOT must span leaves"))
    out=Component[]
    for c in op.components
        aa=Dict(s=>Factor() for s in (:I,:X,:Y,:Z)); bb=Dict(s=>Factor() for s in (:I,:X,:Y,:Z))
        for (p,v) in c.factors[ia]; aa[letter(p,g.a)][p]=v; end
        for (p,v) in c.factors[ib]; bb[letter(p,g.b)][p]=v; end
        for (sa,fa) in aa, (sb,fb) in bb
            (isempty(fa) || isempty(fb)) && continue
            probe=B.word([sa,sb],[g.a,g.b]); mapped,sgn=first(B.expand(g,probe,1.0))
            na=letter(mapped,g.a); nb=letter(mapped,g.b)
            fs=copy(c.factors)
            fs[ia]=Factor(setletter(p,g.a,na)=>v for (p,v) in fa)
            fs[ib]=Factor(setletter(p,g.b,nb)=>v for (p,v) in fb)
            append_checked!(op,out,Component(c.coeff*sgn,fs))
        end
    end
    op.components=out
end
function propagate!(op::Operator,g)
    g isa Union{B.Rotation,B.Clifford} || throw(ArgumentError("unsupported gate $(typeof(g))"))
    if g isa B.Rotation
        B.validword(g.generator,op.n)
        support(g)==0 && return op # identity rotation is a global phase
    else
        max(g.a,g.b)<=op.n || throw(ArgumentError("gate exceeds partition width"))
    end
    ids=affected(op,g)
    isempty(ids) && throw(ArgumentError("gate outside partition"))
    owner=lowest_common_ancestor(op.root,support(g))
    length(ids)>1 && owner.left===nothing && throw(ArgumentError("invalid cross-partition support"))
    if length(ids)==1
        local_gate(op,g,only(ids))
    else
        op.cross_gates+=1
        if g isa B.Rotation
            cross_rotation!(op,g,ids)
        elseif g.kind==:CNOT && length(ids)==2
            cross_cnot!(op,g,ids)
        else
            throw(ArgumentError("unsupported cross-partition gate"))
        end
    end
    op.peak_contributions=max(op.peak_contributions,length(op.components))
    finalize_components!(op)
end
function spatial(c::B.Circuit,terms; leaf_size::Int=2, limits=Limits())
    op=Operator(c.nqubits,leaf_size,terms;limits)
    # Circuit gates are stored in Schrodinger order; Heisenberg propagation reverses them.
    for g in Iterators.reverse(c.gates)
        propagate!(op,g)
    end
    op
end
function sequential(c::B.Circuit,terms; limits=Limits())
    d=Factor()
    for (p,v) in terms; d[p]=get(d,p,0.0)+Float64(v); end
    for g in Iterators.reverse(c.gates)
        out=Factor()
        for (p,v) in d, (q,w) in B.expand(g,p,v)
            w==0 && continue
            out[q]=get(out,q,0.0)+w
            length(out)>limits.max_contributions && throw(ResourceLimit("sequential contribution limit exceeded"))
        end
        d=Factor(p=>v for (p,v) in out if v != 0.0)
    end
    d
end
function expand_global(op::Operator; max_terms::Int=1_000_000)
    out=Factor()
    for c in op.components
        partial=Factor(Word(0,0)=>c.coeff)
        for f in c.factors
            next=Factor()
            length(partial)*length(f)>max_terms && throw(ResourceLimit("test expansion limit exceeded"))
            for (p,v) in partial, (q,w) in f
                r=Word(p.x | q.x,p.z | q.z)
                next[r]=get(next,r,0.0)+v*w
            end
            partial=next
        end
        for (p,v) in partial; out[p]=get(out,p,0.0)+v; end
    end
    Factor(p=>v for (p,v) in out if v != 0.0)
end
function expectation_zero(op::Operator)
    sum((c.coeff*prod((sum(v for (p,v) in f if p.x==0; init=0.0) for f in c.factors); init=1.0) for c in op.components); init=0.0)
end
expectation_zero(d::Factor)=sum(v for (p,v) in d if p.x==0; init=0.0)
stats(op::Operator)=(components=length(op.components),peak_components=op.peak_components,peak_terms=op.peak_terms,peak_contributions=op.peak_contributions,cross_gates=op.cross_gates,before_aggregation=op.before_aggregation,aggregation_ns=op.aggregation_ns)

include("RecursivePauli.jl")

end
