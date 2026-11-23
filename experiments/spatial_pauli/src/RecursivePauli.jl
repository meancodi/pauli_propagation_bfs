# Exact recursive expression DAG. Each internal tensor node owns child expressions;
# sums and gate updates preserve shared, unchanged subtrees.
export RecursiveOperator, recursive_spatial, recursive_expand, recursive_expectation, recursive_stats

abstract type RExpr end
struct RLeaf <: RExpr
    node::Node
    terms::Dict{Word,ComplexF64}
end
struct RTensor <: RExpr
    node::Node
    left::RExpr
    right::RExpr
end
struct RSum <: RExpr
    node::Node
    terms::Vector{RExpr}
end
const MaybeExpr=Union{Nothing,RExpr}
supportnode(x::RExpr)=x.node
function same_expr(a::RExpr,b::RExpr)
    typeof(a)===typeof(b) || return false
    a.node.lo==b.node.lo && a.node.hi==b.node.hi || return false
    a isa RLeaf && return a.terms==b.terms
    a isa RTensor && return same_expr(a.left,b.left) && same_expr(a.right,b.right)
    aa=a::RSum; bb=b::RSum
    length(aa.terms)==length(bb.terms) && all(same_expr(x,y) for (x,y) in zip(aa.terms,bb.terms))
end
function rscale(x::MaybeExpr,c::Number)
    (x===nothing || c==0) && return nothing
    c==1 && return x
    if x isa RLeaf
        return RLeaf(x.node,Dict(p=>ComplexF64(c*v) for (p,v) in x.terms if c*v!=0))
    elseif x isa RTensor
        return RTensor(x.node,rscale(x.left,c)::RExpr,x.right)
    else
        return rsum(x.node,[rscale(t,c) for t in x.terms])
    end
end
function rtensor(node::Node,left::MaybeExpr,right::MaybeExpr)
    (left===nothing || right===nothing) && return nothing
    RTensor(node,left,right)
end
function rsum(node::Node,terms)
    flat=RExpr[]
    for t in terms
        t===nothing && continue
        if t isa RSum
            append!(flat,t.terms)
        else
            push!(flat,t)
        end
    end
    isempty(flat) && return nothing
    # Common child factoring is an exact internal-node aggregation.
    changed=true
    while changed
        changed=false
        for i in eachindex(flat)
            i>=length(flat) && break
            for j in i+1:length(flat)
                a=flat[i]; b=flat[j]
                if a isa RLeaf && b isa RLeaf
                    d=copy(a.terms)
                    for (p,v) in b.terms; d[p]=get(d,p,0im)+v; end
                    flat[i]=RLeaf(node,Dict(p=>v for (p,v) in d if v!=0))
                elseif a isa RTensor && b isa RTensor && same_expr(a.right,b.right)
                    left=rsum(a.left.node,[a.left,b.left])
                    flat[i]=rtensor(node,left,a.right)
                elseif a isa RTensor && b isa RTensor && same_expr(a.left,b.left)
                    right=rsum(a.right.node,[a.right,b.right])
                    flat[i]=rtensor(node,a.left,right)
                else
                    continue
                end
                deleteat!(flat,j); changed=true; break
            end
            changed && break
        end
    end
    filter!(t->!(t isa RLeaf && isempty(t.terms)),flat)
    isempty(flat) ? nothing : length(flat)==1 ? only(flat) : RSum(node,flat)
end
function rword(node::Node,p::Word,c::ComplexF64)
    if node.left===nothing
        return RLeaf(node,Dict(Word(p.x & mask64(node),p.z & mask64(node))=>c))
    end
    RTensor(node,rword(node.left,p,c),rword(node.right,p,1.0+0im))
end
mutable struct RecursiveOperator
    n::Int
    root::Node
    expr::MaybeExpr
    reuse_events::Int
    internal_aggregations::Int
    cross_by_level::Dict{Int,Int}
    peak_nodes::Int
    peak_terms::Int
    limits::Limits
end
function RecursiveOperator(n::Int,leaf_size::Int,terms;limits=Limits())
    root=partition(n,leaf_size)
    exprs=RExpr[]
    for (p,c) in terms
        B.validword(p,n); isfinite(c) || throw(ArgumentError("nonfinite coefficient"))
        c==0 && continue
        push!(exprs,rword(root,p,ComplexF64(c)))
    end
    expr=rsum(root,exprs)
    op=RecursiveOperator(n,root,expr,0,0,Dict{Int,Int}(),0,0,limits)
    check_recursive!(op)
end
function rsplit(x::RLeaf,gen::Word)
    localgen=Word(gen.x & mask64(x.node),gen.z & mask64(x.node))
    a=Dict{Word,ComplexF64}(); b=Dict{Word,ComplexF64}()
    for (p,v) in x.terms
        (B.commutes(localgen,p) ? a : b)[p]=v
    end
    (isempty(a) ? nothing : RLeaf(x.node,a),isempty(b) ? nothing : RLeaf(x.node,b))
end
function rsplit(x::RTensor,gen::Word)
    lc,la=rsplit(x.left,gen); rc,ra=rsplit(x.right,gen)
    comm=rsum(x.node,[rtensor(x.node,lc,rc),rtensor(x.node,la,ra)])
    anti=rsum(x.node,[rtensor(x.node,lc,ra),rtensor(x.node,la,rc)])
    comm,anti
end
function rsplit(x::RSum,gen::Word)
    pieces=[rsplit(t,gen) for t in x.terms]
    rsum(x.node,[p[1] for p in pieces]),rsum(x.node,[p[2] for p in pieces])
end
function rmultiply(x::RLeaf,gen::Word)
    localgen=Word(gen.x & mask64(x.node),gen.z & mask64(x.node))
    d=Dict{Word,ComplexF64}()
    for (p,v) in x.terms
        q,phase=B.product(localgen,p)
        d[q]=get(d,q,0im)+v*phase
    end
    RLeaf(x.node,d)
end
rmultiply(x::RTensor,gen::Word)=RTensor(x.node,rmultiply(x.left,gen),rmultiply(x.right,gen))
rmultiply(x::RSum,gen::Word)=rsum(x.node,[rmultiply(t,gen) for t in x.terms])
function rproject(x::RLeaf,q::Int,s::Symbol)
    d=Dict(p=>v for (p,v) in x.terms if letter(p,q)==s)
    isempty(d) ? nothing : RLeaf(x.node,d)
end
function rproject(x::RTensor,q::Int,s::Symbol)
    if q<=x.left.node.hi
        return rtensor(x.node,rproject(x.left,q,s),x.right)
    end
    rtensor(x.node,x.left,rproject(x.right,q,s))
end
rproject(x::RSum,q::Int,s::Symbol)=rsum(x.node,[rproject(t,q,s) for t in x.terms])
function rsetletter(x::RLeaf,q::Int,s::Symbol)
    RLeaf(x.node,Dict(setletter(p,q,s)=>v for (p,v) in x.terms))
end
function rsetletter(x::RTensor,q::Int,s::Symbol)
    q<=x.left.node.hi ? RTensor(x.node,rsetletter(x.left,q,s),x.right) : RTensor(x.node,x.left,rsetletter(x.right,q,s))
end
rsetletter(x::RSum,q::Int,s::Symbol)=rsum(x.node,[rsetletter(t,q,s) for t in x.terms])
function rlocal(x::RLeaf,g)
    d=Dict{Word,ComplexF64}()
    for (p,v) in x.terms, (q,w) in B.expand(g,p,1.0)
        w==0 && continue
        d[q]=get(d,q,0im)+v*w
    end
    RLeaf(x.node,Dict(p=>v for (p,v) in d if v!=0))
end
function rcross_cnot(x::RTensor,g::B.Clifford)
    pieces=RExpr[]
    for sa in (:I,:X,:Y,:Z), sb in (:I,:X,:Y,:Z)
        left=rproject(x.left,g.a<=x.left.node.hi ? g.a : g.b,g.a<=x.left.node.hi ? sa : sb)
        right=rproject(x.right,g.a>x.left.node.hi ? g.a : g.b,g.a>x.left.node.hi ? sa : sb)
        (left===nothing || right===nothing) && continue
        probe=B.word([sa,sb],[g.a,g.b]); mapped,sgn=first(B.expand(g,probe,1.0))
        if g.a<=x.left.node.hi
            left=rsetletter(left,g.a,letter(mapped,g.a)); right=rsetletter(right,g.b,letter(mapped,g.b))
        else
            left=rsetletter(left,g.b,letter(mapped,g.b)); right=rsetletter(right,g.a,letter(mapped,g.a))
        end
        push!(pieces,rscale(RTensor(x.node,left,right),sgn)::RExpr)
    end
    rsum(x.node,pieces)
end
function rapply(x::RLeaf,g,op::RecursiveOperator)
    rlocal(x,g)
end
function rapply(x::RSum,g,op::RecursiveOperator)
    op.internal_aggregations+=1
    rsum(x.node,[rapply(t,g,op) for t in x.terms])
end
function rapply(x::RTensor,g,op::RecursiveOperator)
    s=support(g)
    if (s & ~mask64(x.left.node))==0
        op.reuse_events+=1
        return rtensor(x.node,rapply(x.left,g,op),x.right)
    elseif (s & ~mask64(x.right.node))==0
        op.reuse_events+=1
        return rtensor(x.node,x.left,rapply(x.right,g,op))
    end
    if g isa B.Rotation
        comm,anti=rsplit(x,g.generator)
        return rsum(x.node,[comm,rscale(anti,g.cosine),anti===nothing ? nothing : rscale(rmultiply(anti,g.generator),im*g.sine)])
    end
    g.kind==:CNOT || throw(ArgumentError("unsupported cross-tree gate"))
    rcross_cnot(x,g)
end
count_nodes_to_leaf(n::Node)=n.left===nothing ? 0 : 1+max(count_nodes_to_leaf(n.left),count_nodes_to_leaf(n.right))
function exprsize(x::MaybeExpr)
    x===nothing && return (0,0)
    x isa RLeaf && return (1,length(x.terms))
    kids=x isa RTensor ? (x.left,x.right) : x.terms
    counts=map(exprsize,kids)
    (1+sum(y[1] for y in counts),sum(y[2] for y in counts))
end
function check_recursive!(op::RecursiveOperator)
    nodes,terms=exprsize(op.expr)
    op.peak_nodes=max(op.peak_nodes,nodes); op.peak_terms=max(op.peak_terms,terms)
    nodes>op.limits.max_contributions && throw(ResourceLimit("recursive expression node limit exceeded"))
    terms>op.limits.max_local_terms && throw(ResourceLimit("recursive local term limit exceeded"))
    op.limits.max_representation_bytes!=typemax(Int) && Base.summarysize(op.expr)>op.limits.max_representation_bytes && throw(ResourceLimit("recursive byte estimate exceeded"))
    op
end
function recursive_spatial(c::B.Circuit,terms;leaf_size::Int=2,limits=Limits())
    op=RecursiveOperator(c.nqubits,leaf_size,terms;limits)
    for g in Iterators.reverse(c.gates)
        op.expr===nothing && break
        owner=lowest_common_ancestor(op.root,support(g))
        if owner.left!==nothing
            level=count_nodes_to_leaf(owner)
            op.cross_by_level[level]=get(op.cross_by_level,level,0)+1
        end
        op.expr=rapply(op.expr,g,op)
        check_recursive!(op)
    end
    op
end
function recursive_expand(x::RLeaf;max_terms::Int=1_000_000)
    copy(x.terms)
end
function recursive_expand(x::RTensor;max_terms::Int=1_000_000)
    a=recursive_expand(x.left;max_terms); b=recursive_expand(x.right;max_terms)
    length(a)*length(b)>max_terms && throw(ResourceLimit("recursive test expansion limit exceeded"))
    d=Dict{Word,ComplexF64}()
    for (p,v) in a, (q,w) in b
        r=Word(p.x|q.x,p.z|q.z)
        d[r]=get(d,r,0im)+v*w
    end
    d
end
function recursive_expand(x::RSum;max_terms::Int=1_000_000)
    d=Dict{Word,ComplexF64}()
    for t in x.terms, (p,v) in recursive_expand(t;max_terms)
        d[p]=get(d,p,0im)+v
        length(d)>max_terms && throw(ResourceLimit("recursive test expansion limit exceeded"))
    end
    d
end
recursive_expand(op::RecursiveOperator;max_terms::Int=1_000_000)=op.expr===nothing ? Dict{Word,ComplexF64}() : recursive_expand(op.expr;max_terms)
recursive_expectation(x::RLeaf)=sum((v for (p,v) in x.terms if p.x==0);init=0im)
recursive_expectation(x::RTensor)=recursive_expectation(x.left)*recursive_expectation(x.right)
recursive_expectation(x::RSum)=sum((recursive_expectation(t) for t in x.terms);init=0im)
recursive_expectation(op::RecursiveOperator)=op.expr===nothing ? 0.0 : real(recursive_expectation(op.expr))
function per_node_components(x::MaybeExpr)
    counts=Dict{Tuple{Int,Int},Int}()
    function visit(y::MaybeExpr)
        y===nothing && return
        k=(y.node.lo,y.node.hi)
        counts[k]=max(get(counts,k,0),y isa RSum ? length(y.terms) : 1)
        if y isa RTensor
            visit(y.left); visit(y.right)
        elseif y isa RSum
            for t in y.terms; visit(t); end
        end
    end
    visit(x)
    counts
end
recursive_stats(op::RecursiveOperator)=(nodes=exprsize(op.expr)[1],peak_nodes=op.peak_nodes,peak_terms=op.peak_terms,reuse_events=op.reuse_events,internal_aggregations=op.internal_aggregations,cross_by_level=copy(op.cross_by_level),per_node_components=per_node_components(op.expr))
