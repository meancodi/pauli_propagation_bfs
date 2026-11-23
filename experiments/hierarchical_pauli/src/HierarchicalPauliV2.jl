module HierarchicalPauliV2

using Random, LinearAlgebra, Statistics
const B = Main.BestFirstPauli
const PW = B.PauliWord
const Terms = Dict{PW,Float64}

export V2Config, V2Metrics, V2ColumnMap, V2Node, MaterializedNode,
       SequentialNode, OnDemandNode, build_plan, apply_plan, v2_propagate,
       local_words, truncate_column, node_support, node_entries, plan_summary,
       map_metrics, reset_metrics, make_config

bitword(a,q)=!iszero(a & (UInt64(1)<<(q-1)))
function merge_terms(terms)
    d=Terms()
    for (p,c) in terms
        c==0 && continue
        d[p]=get(d,p,0.0)+c
    end
    filter!(x->x.second != 0.0,d); d
end
function aggregate_terms(terms; tol=1e-12)
    d=Terms(); duplicates=0; cancellations=0
    for (p,c) in terms
        if haskey(d,p); duplicates+=1; v=d[p]+c; abs(v)<=tol && (cancellations+=1); d[p]=v
        else; d[p]=c end
    end
    filter!(x->x.second != 0.0,d); d,duplicates,cancellations
end
function local_words(support)
    out=PW[]
    for code in 0:(4^length(support)-1)
        x=z=UInt64(0); t=code
        for q in support
            v=mod(t,4); t ÷= 4; b=UInt64(1)<<(q-1)
            v==1 && (x|=b); v==2 && (z|=b); v==3 && ((x|=b);(z|=b))
        end
        push!(out,PW(x,z))
    end
    out
end

Base.@kwdef mutable struct V2Config
    truncation_mode::Symbol=:none
    coefficient_threshold::Float64=0.0
    max_column_terms::Int=typemax(Int)
    selected_levels::Set{Int}=Set{Int}()
    size_trigger_contributions::Int=typemax(Int)
    size_trigger_entries::Int=typemax(Int)
    max_support::Int=typemax(Int)
    max_columns::Int=typemax(Int)
    max_entries::Int=typemax(Int)
    max_temp_contributions::Int=typemax(Int)
    max_live_terms::Int=typemax(Int)
    max_seconds::Float64=Inf
end
make_config(;kwargs...)=V2Config(;kwargs...)

Base.@kwdef mutable struct V2Metrics
    total_generated::Int=0
    current_generated::Int=0
    duplicate_contributions::Int=0
    true_cancellations::Int=0
    below_threshold::Dict{Float64,Int}=Dict{Float64,Int}()
    discarded_mass::Float64=0.0
    discarded_terms::Int=0
    map_columns::Int=0
    map_entries::Int=0
    peak_temp_contributions::Int=0
    peak_map_entries::Int=0
    leaf_seconds::Float64=0.0
    composition_seconds::Float64=0.0
    truncation_seconds::Float64=0.0
    application_seconds::Float64=0.0
    tree_seconds::Float64=0.0
    internal_truncations::Int=0
    materialized_nodes::Int=0
    sequential_nodes::Int=0
    fallback_nodes::Int=0
    recursion_nodes::Int=0
    resource_events::Int=0
    cache_hits::Int=0
    cache_misses::Int=0
    unavailable_columns::Int=0
end
reset_metrics()=V2Metrics()

struct V2ColumnMap
    columns::Dict{PW,Terms}
    errors::Dict{PW,Float64}
    support::Vector{Int}
    level::Int
    metrics::V2Metrics
end
abstract type V2Node end
struct MaterializedNode <: V2Node
    map::V2ColumnMap
    gates::Vector{B.SupportedGate}
end
struct SequentialNode <: V2Node
    left::V2Node
    right::V2Node
    support::Vector{Int}
    reason::Symbol
end
mutable struct OnDemandNode <: V2Node
    gates::Vector{B.SupportedGate}
    support::Vector{Int}
    cache::Dict{PW,Tuple{Terms,Float64}}
    metrics::V2Metrics
end

node_support(n::MaterializedNode)=n.map.support
node_support(n::OnDemandNode)=n.support
node_support(n::SequentialNode)=n.support
node_entries(n::MaterializedNode)=sum(length(v) for v in values(n.map.columns))
node_entries(::OnDemandNode)=0
node_entries(n::SequentialNode)=node_entries(n.left)+node_entries(n.right)
function map_metrics(n::V2Node)
    n isa MaterializedNode && return n.map.metrics
    n isa OnDemandNode && return n.metrics
    m=V2Metrics(); merge_metrics!(m,map_metrics(n.left)); merge_metrics!(m,map_metrics(n.right)); m
end
function merge_metrics!(a::V2Metrics,b::V2Metrics)
    for f in fieldnames(V2Metrics)
        v=getfield(a,f); w=getfield(b,f)
        if v isa Dict
            for (k,x) in w; v[k]=get(v,k,0)+x end
        elseif v isa Int; setfield!(a,f,v+w)
        elseif v isa Float64; setfield!(a,f,v+w)
        end
    end
    a
end

function gate_map(g,inputs,m::V2Metrics)
    out=Tuple{PW,Float64}[]
    for (p,c) in inputs
        for (q,cq) in B.expand(g,p,c); push!(out,(q,cq)); end
    end
    m.total_generated+=length(out); m.current_generated=length(out); m.peak_temp_contributions=max(m.peak_temp_contributions,length(out))
    d,du,ca=aggregate_terms(out); m.duplicate_contributions+=du; m.true_cancellations+=ca
    out,d
end
function apply_gates(gates,p,m::V2Metrics)
    d=Terms(p=>1.0)
    for g in reverse(gates); _,d=gate_map(g,d,m); end
    d
end
function support_of(gates)
    qs=Int[]
    for g in gates
        if g isa B.Rotation
            for q in 1:64; bitword(g.generator.x|g.generator.z,q)&&push!(qs,q) end
        else
            push!(qs,g.a); g.kind==:CNOT&&push!(qs,g.b)
        end
    end
    sort(unique(qs))
end
function apply_local_map(map::V2ColumnMap,p,c=1.0)
    mask=zero(UInt64); for q in map.support; mask|=UInt64(1)<<(q-1) end
    localp=PW(p.x&mask,p.z&mask); col=get(map.columns,localp,nothing)
    col===nothing && throw(ArgumentError("unavailable column for materialized map"))
    Terms(PW((q.x&mask)|(p.x&~mask),(q.z&mask)|(p.z&~mask))=>c*v for (q,v) in col)
end
function choose_truncation(config,level,isroot,generated,entries)
    mode=config.truncation_mode
    mode==:none && return false
    mode==:root_only && return isroot
    mode==:all_internal && return true
    mode==:alternate_levels && return iseven(level)
    mode==:selected_levels && return level in config.selected_levels
    mode==:size_triggered && return generated>config.size_trigger_contributions || entries>config.size_trigger_entries
    throw(ArgumentError("unknown truncation mode $mode"))
end
function truncate_column(d::Terms,err::Float64,config::V2Config,m::V2Metrics)
    t=time_ns(); out=copy(d); extra=0.0; nd=0
    for (p,c) in sort(collect(out);by=x->(abs(x[2]),x[1].x,x[1].z))
        if abs(c)<config.coefficient_threshold; extra+=abs(c); delete!(out,p); nd+=1 end
    end
    if length(out)>config.max_column_terms
        keep=sort(collect(out);by=x->(-abs(x[2]),x[1].x,x[1].z))
        for (p,c) in keep[config.max_column_terms+1:end]; extra+=abs(c); delete!(out,p); nd+=1 end
    end
    if config.coefficient_threshold>0
        m.below_threshold[config.coefficient_threshold]=get(m.below_threshold,config.coefficient_threshold,0)+nd
    end
    m.discarded_mass+=extra; m.discarded_terms+=nd; m.internal_truncations+=1
    m.truncation_seconds+=(time_ns()-t)/1e9
    out,err+extra
end
function build_leaf(gates,config,m)
    support=support_of(gates); length(support)<=config.max_support || return nothing
    nc=try 4^length(support) catch; typemax(Int) end
    nc<=config.max_columns || return nothing
    t=time_ns(); cols=Dict{PW,Terms}(); errs=Dict{PW,Float64}()
    for p in local_words(support)
        d=apply_gates(gates,p,m); cols[p]=d; errs[p]=0.0
        m.peak_map_entries=max(m.peak_map_entries,sum(length(v) for v in values(cols)))
    end
    entries=sum((length(v) for v in values(cols));init=0)
    entries>config.max_entries && return nothing
    m.leaf_seconds+=(time_ns()-t)/1e9; m.map_columns+=length(cols); m.map_entries+=entries; m.materialized_nodes+=1
    V2ColumnMap(cols,errs,support,0,m)
end
function compose_maps(left::V2ColumnMap,right::V2ColumnMap,config,level,isroot,m)
    support=sort(unique(vcat(left.support,right.support))); nc=try 4^length(support) catch; typemax(Int) end
    nc<=config.max_columns || return nothing
    t=time_ns(); cols=Dict{PW,Terms}(); errs=Dict{PW,Float64}(); temp=0
    for p in local_words(support)
        rcol=apply_local_map(right,p); terms=Tuple{PW,Float64}[]; eright=0.0
        rlocal=PW[]
        for (mid,cr) in rcol
            push!(rlocal,mid); eright+=abs(cr)*get(right.errors,PW(mid.x & sum(UInt64(1)<<(q-1) for q in right.support),mid.z & sum(UInt64(1)<<(q-1) for q in right.support)),0.0)
            for (q,c) in apply_local_map(left,mid); push!(terms,(q,cr*c)) end
        end
        temp=max(temp,length(terms)); m.total_generated+=length(terms); m.current_generated=length(terms); m.peak_temp_contributions=max(m.peak_temp_contributions,temp)
        d,du,ca=aggregate_terms(terms); m.duplicate_contributions+=du; m.true_cancellations+=ca
        e=eright
        for mid in rlocal
            mask=zero(UInt64); for q in left.support; mask|=UInt64(1)<<(q-1) end
            lm=PW(mid.x&mask,mid.z&mask)
            # Coefficients are recomputed from the right column below; this term
            # is intentionally conservative and uses the column 1-norm.
            e+=sum(abs(c)*get(left.errors,q,0.0) for (q,c) in apply_local_map(left,mid))
        end
        if choose_truncation(config,level,isroot,length(terms),sum((length(v) for v in values(cols));init=0))
            d,e=truncate_column(d,e,config,m)
        end
        cols[p]=d; errs[p]=e
    end
    newentries=sum((length(v) for v in values(cols));init=0)
    if newentries>config.max_entries
        m.resource_events+=1; return nothing
    end
    m.composition_seconds+=(time_ns()-t)/1e9; m.map_columns+=length(cols); m.map_entries+=newentries; m.peak_map_entries=max(m.peak_map_entries,newentries)
    V2ColumnMap(cols,errs,support,level,m)
end

function estimated_support(gates,i,j); support_of(gates[i:j]) end
function choose_split(gates,strategy,config)
    n=length(gates); n<=1 && return 0
    strategy==:left_deep && return n-1
    strategy==:balanced && return n÷2
    candidates=1:n-1; mid=n/2
    scores=map(candidates) do k
        s=length(union(estimated_support(gates,1,k),estimated_support(gates,k+1,n)))
        potential=try 4^s catch; typemax(Int) end
        (s>config.max_support || potential>config.max_columns ? 1 : 0, max(s,potential>config.max_columns ? typemax(Int) : potential),abs(k-mid))
    end
    candidates[argmin(scores)]
end
function build_plan(gates;strategy=:balanced,config=V2Config(),level=0,isroot=true,m=V2Metrics())
    m.recursion_nodes+=1
    if length(gates)==1
        leaf=build_leaf(gates,config,m)
        leaf===nothing && (m.resource_events+=1; m.fallback_nodes+=1; return OnDemandNode(copy(gates),support_of(gates),Dict{PW,Tuple{Terms,Float64}}(),m))
        return MaterializedNode(leaf,copy(gates))
    end
    k=choose_split(gates,strategy,config); left=build_plan(gates[1:k];strategy=config.max_support<typemax(Int) && strategy==:support_constrained ? :support_constrained : strategy,config=config,level=level+1,isroot=false,m=m); right=build_plan(gates[k+1:end];strategy=strategy,config=config,level=level+1,isroot=false,m=m)
    if left isa MaterializedNode && right isa MaterializedNode
        parent=compose_maps(left.map,right.map,config,level+1,isroot,m)
        parent!==nothing && return MaterializedNode(parent,vcat(left.gates,right.gates))
    end
    m.sequential_nodes+=1; m.fallback_nodes+=1; m.resource_events+=1
    SequentialNode(left,right,sort(unique(vcat(node_support(left),node_support(right)))), :resource_limit)
end

function apply_node(n::MaterializedNode,d,err)
    out=Terms(); e=err
    for (p,c) in d
        col=apply_local_map(n.map,p,c); for (q,v) in col; out[q]=get(out,q,0.0)+v end
        localp=PW(p.x & sum(UInt64(1)<<(q-1) for q in n.map.support),p.z & sum(UInt64(1)<<(q-1) for q in n.map.support)); e+=abs(c)*get(n.map.errors,localp,0.0)
    end
    merge_terms(out),e
end
function apply_node(n::OnDemandNode,d,err)
    out=Terms(); e=err; mask=zero(UInt64); for q in n.support; mask|=UInt64(1)<<(q-1) end
    for (p,c) in d
        lp=PW(p.x&mask,p.z&mask); entry=get(n.cache,lp,nothing)
        if entry===nothing
            n.metrics.cache_misses+=1; col=apply_gates(n.gates,lp,n.metrics); n.cache[lp]=(col,0.0)
        else
            n.metrics.cache_hits+=1; col=entry[1]
        end
        for (q,v) in col; out[PW((q.x&mask)|(p.x&~mask),(q.z&mask)|(p.z&~mask))]=get(out,PW((q.x&mask)|(p.x&~mask),(q.z&mask)|(p.z&~mask)),0.0)+c*v end
    end
    merge_terms(out),e
end
function apply_node(n::SequentialNode,d,err)
    d1,e1=apply_node(n.right,d,err); apply_node(n.left,d1,e1)
end
function apply_plan(n,d::Terms=Terms(),err=0.0)
    isempty(d) && return (Terms(),err)
    t=time_ns(); out,e=apply_node(n,d,err); map_metrics(n).application_seconds+=(time_ns()-t)/1e9; out,e
end
function v2_propagate(gates,obs;strategy=:balanced,config=V2Config())
    started=time_ns(); m=V2Metrics(); plan=build_plan(gates;strategy=strategy,config=config,m=m,isroot=true); terms,eb=apply_plan(plan,Terms(obs=>1.0),0.0)
    m.tree_seconds=(time_ns()-started)/1e9
    (;plan,terms,expectation=sum((c for (p,c) in terms if p.x==0);init=0.0),observed_error=missing,certified_bound=eb,
      metrics=m,status=m.resource_events>0 ? :resource_limited : :complete,reason=m.resource_events>0 ? :fallback : :finished)
end
function plan_summary(n::V2Node)
    n isa MaterializedNode && return (kind=:materialized,support=length(n.map.support),entries=node_entries(n))
    n isa OnDemandNode && return (kind=:on_demand,support=length(n.support),entries=length(n.cache))
    (kind=:sequential,support=length(n.support),left=plan_summary(n.left),right=plan_summary(n.right),reason=n.reason)
end
end
