using Dates, SHA, Statistics
include(joinpath(@__DIR__,"run_benchmarks.jl"))
const PP=B.PP
const V2ROOT=normpath(joinpath(@__DIR__,"..","results","v2"))
const V2LIMITS=Limits(max_components=4_000,max_local_terms=60_000,max_contributions=12_000,max_aggregation_passes=4)
const REFLIMITS=Limits(max_components=100_000,max_local_terms=100_000,max_contributions=100_000)
const FAMILY_CODE=Dict(:pairs=>11,:ising=>13,:periodic_ising=>17,:weak=>19,:cross=>23,:clifford=>29)

function v2circuit(n,family,layers,seed,coupling)
    periodic=family==:periodic_ising
    circuit(n,periodic ? :ising : family,layers,seed;coupling,periodic)
end
function v2observable(n,kind)
    kind==:factorized && return [B.word([a,b],[1,n])=>1.0 for a in (:X,:Z), b in (:X,:Z)]
    observable(n,kind)
end
function fingerprint(c,terms)
    parts=String[]
    for g in c.gates
        if g isa B.Rotation
            push!(parts,"R:$(g.generator.x):$(g.generator.z):$(repr(g.theta))")
        else
            push!(parts,"C:$(g.kind):$(g.a):$(g.b)")
        end
    end
    for (p,v) in terms; push!(parts,"O:$(p.x):$(p.z):$(repr(v))"); end
    bytes2hex(sha256(join(parts,"|")))
end
function upstream_exact(c,terms)
    gates,theta,_=B.upstream_inputs(c,first(terms).first)
    psum=PP.PauliSum(c.nqubits)
    for (p,v) in terms
        PP.add!(psum,PP.PauliString(c.nqubits,B.symbols(p,c.nqubits),collect(1:c.nqubits),v))
    end
    PP.propagate(gates,psum,theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
end
function v2run(method,c,terms;profile=false)
    t0=time_ns(); init_ns=0; gate_ns=0; agg_ns=0; peak=0; obj=nothing
    if method=="S0"
        d=SpatialPauli.Factor()
        for (p,v) in terms; d[p]=get(d,p,0.0)+v; end
        init_ns=time_ns()-t0; profile && (peak=Base.summarysize(d))
        t=time_ns()
        for g in Iterators.reverse(c.gates)
            next=SpatialPauli.Factor()
            for (p,v) in d, (q,w) in B.expand(g,p,v)
                w==0 && continue
                next[q]=get(next,q,0.0)+w
                length(next)>REFLIMITS.max_contributions && throw(ResourceLimit("S0 contribution limit"))
            end
            d=SpatialPauli.Factor(p=>v for (p,v) in next if v!=0)
            profile && (peak=max(peak,Base.summarysize(d)))
        end
        gate_ns=time_ns()-t; obj=d
    elseif method=="S0_optimized"
        obj=upstream_exact(c,terms)
        init_ns=0; gate_ns=time_ns()-t0
        profile && (peak=Base.summarysize(obj))
    elseif method in ("S1","S2")
        leaf=method=="S1" ? cld(c.nqubits,2) : 2
        obj=Operator(c.nqubits,leaf,terms;limits=V2LIMITS)
        init_ns=time_ns()-t0; agg_before=obj.aggregation_ns
        profile && (peak=Base.summarysize(obj))
        t=time_ns()
        for g in Iterators.reverse(c.gates)
            SpatialPauli.propagate!(obj,g)
            profile && (peak=max(peak,Base.summarysize(obj)))
        end
        gate_ns=time_ns()-t; agg_ns=obj.aggregation_ns-agg_before
    elseif method=="S3"
        obj=RecursiveOperator(c.nqubits,2,terms;limits=V2LIMITS)
        init_ns=time_ns()-t0
        profile && (peak=Base.summarysize(obj))
        t=time_ns()
        for g in Iterators.reverse(c.gates)
            obj.expr===nothing && break
            obj.expr=SpatialPauli.rapply(obj.expr,g,obj)
            SpatialPauli.check_recursive!(obj)
            profile && (peak=max(peak,Base.summarysize(obj)))
        end
        gate_ns=time_ns()-t
    else
        throw(ArgumentError("unknown method $method"))
    end
    total_ns=time_ns()-t0
    expectation=method=="S0_optimized" ? PP.overlapwithzero(obj) : method=="S3" ? recursive_expectation(obj) : expectation_zero(obj)
    components=method in ("S1","S2") ? length(obj.components) : 0
    peak_components=method in ("S1","S2") ? obj.peak_components : 0
    peak_terms=method in ("S1","S2") ? obj.peak_terms : method=="S3" ? obj.peak_terms : 0
    peak_contributions=method in ("S1","S2") ? obj.peak_contributions : method=="S3" ? obj.peak_nodes : 0
    reuse=method=="S3" ? obj.reuse_events : 0
    internal=method=="S3" ? obj.internal_aggregations : 0
    cross=method=="S3" ? repr(obj.cross_by_level) : method in ("S1","S2") ? string(obj.cross_gates) : ""
    per_node=method=="S3" ? repr(recursive_stats(obj).per_node_components) : ""
    (;obj,init_ns,gate_ns,agg_ns,total_ns,peak,expectation,components,peak_components,peak_terms,peak_contributions,reuse,internal,cross,per_node)
end
function v2cases(mode)
    if mode=="pilot"
        return [(n,f,1,k,1) for n in (4,8,12), f in (:pairs,:cross), k in (:single,:factorized)]
    elseif mode=="final"
        return vcat(
            vec([(8,f,1,k,1) for f in (:pairs,:ising,:periodic_ising,:weak,:cross,:clifford), k in (:zz,:factorized)]),
            [(8,f,2,:zz,1) for f in (:pairs,:ising,:weak,:clifford)],
            vec([(12,f,1,k,1) for f in (:pairs,:ising,:periodic_ising,:weak,:clifford), k in (:single,:factorized)]),
            [(16,f,1,:factorized,1) for f in (:pairs,:clifford)])
    elseif mode=="targeted"
        return [(8,f,d,k,s) for f in (:pairs,:cross), d in (1,2), k in (:zz,:factorized), s in 1:10]
    elseif mode=="seeded"
        return vcat([(8,:pairs,1,:factorized,s) for s in 1:10],
                    [(8,:cross,1,:factorized,s) for s in 1:10])
    elseif mode=="coupling"
        return [(8,:weak,1,:factorized,j) for j in 1:4]
    elseif mode=="diagnostic"
        return [(n,f,1,:factorized,1) for n in (4,8), f in (:pairs,:cross)]
    end
    throw(ArgumentError("unknown campaign mode"))
end
function writecsv(io,fields)
    println(io,join(("\""*replace(string(x),'"'=>"''")*"\"" for x in fields),","))
end
function main_v2()
    mode=get(ENV,"SPATIAL_V2_MODE","pilot")
    root=joinpath(V2ROOT,"raw"); mkpath(root)
    path=joinpath(root,"v2_"*mode*"_"*Dates.format(now(),"yyyymmdd_HHMMSS")*".csv")
    header=("version","family","n","layers","gate_count","observable","coupling","periodic","seed","fingerprint","partition_scheme","leaf_size","method","status","median_total_s","min_total_s","max_total_s","median_init_s","median_gate_s","median_gate_excluding_aggregation_s","median_aggregation_s","peak_retained_bytes","allocated_bytes","components","peak_components","peak_local_terms","peak_intermediate","subtree_reuse","internal_aggregations","cross_levels","per_node_components","expectation","expectation_error","julia","os","cpu","threads","message")
    open(path,"w") do io; writecsv(io,header); end
    println("Writing ",path)
    cases=v2cases(mode)
    maxcases=parse(Int,get(ENV,"SPATIAL_V2_MAX_CASES",string(length(cases))))
    methods=("S0","S0_optimized","S1","S2","S3")
    for (caseidx,(n,family,layers,kind,seedidx)) in enumerate(cases)
        caseidx>maxcases && break
        coupling=mode=="coupling" ? seedidx-1 : (family in (:weak,:cross) ? layers : 0)
        seed=100_000+FAMILY_CODE[family]*10_000+n*100+layers*10+(mode=="coupling" ? 1 : seedidx)
        c=v2circuit(n,family,layers,seed,coupling); terms=v2observable(n,kind)
        fp=fingerprint(c,terms); periodic=family==:periodic_ising
        ref=nothing; ref_error=""
        try
            ref=sequential(c,terms;limits=REFLIMITS)
        catch e
            ref_error=sprint(showerror,e)
        end
        for method in methods
            status="completed"; msg=""; samples=Any[]; bytes=0; memory=nothing
            if ref===nothing
                status="skipped"; msg="reference unavailable: "*ref_error
            else
                try
                    v2run(method,c,terms)
                    for _ in 1:5
                        GC.gc()
                        b=@allocated r=v2run(method,c,terms)
                        bytes=max(bytes,b); push!(samples,r)
                    end
                    memory=v2run(method,c,terms;profile=true)
                catch e
                    status=e isa ResourceLimit ? "resource_limited" : "error"
                    msg=sprint(showerror,e)
                end
            end
            leaf=method=="S1" ? cld(n,2) : method in ("S2","S3") ? 2 : 0
            scheme=method=="S3" ? "recursive_balanced" : method=="S2" ? "flat_aligned_leaves" : method=="S1" ? "flat_bipartition" : "global"
            row=Any["2.2",family,n,layers,length(c.gates),kind,coupling,periodic,seed,fp,scheme,leaf,method,status]
            if status=="completed"
                times=sort([r.total_ns for r in samples]); med=times[3]
                r=samples[findfirst(x->x.total_ns==med,samples)]
                err=abs(r.expectation-expectation_zero(ref))
                if err>1e-8
                    status="error"; row[end]=status; msg="expectation mismatch $err"
                end
                agg=method in ("S1","S2") ? r.agg_ns/1e9 : ""
                excl=method in ("S1","S2") ? (r.gate_ns-r.agg_ns)/1e9 : ""
                append!(row,[med/1e9,times[1]/1e9,times[end]/1e9,r.init_ns/1e9,r.gate_ns/1e9,excl,agg,memory.peak,bytes,r.components,r.peak_components,r.peak_terms,r.peak_contributions,r.reuse,r.internal,r.cross,r.per_node,r.expectation,err])
            else
                append!(row,fill("",19))
            end
            append!(row,[VERSION,Sys.KERNEL,Sys.CPU_NAME,Threads.nthreads(),msg])
            open(path,"a") do io; writecsv(io,row); end
        end
    end
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_v2()
end
