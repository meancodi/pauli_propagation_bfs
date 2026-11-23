using Random, Statistics, Dates
include(joinpath(@__DIR__,"..","src","SpatialPauli.jl"))
using .SpatialPauli
const B=SpatialPauli.B
const PW=B.PauliWord
const lim=Limits(max_components=5_000,max_local_terms=100_000,max_contributions=20_000)

function circuit(n,family,depth,seed; coupling=1,periodic=false)
    rng=MersenneTwister(seed); gates=B.SupportedGate[]
    for layer in 1:depth
        if family==:pairs || family==:weak || family==:cross
            for q in 1:n
                push!(gates,B.Rotation(B.word([:X],[q]),0.23+0.08rand(rng)))
                push!(gates,B.Rotation(B.word([:Z],[q]),0.14+0.07rand(rng)))
            end
            for q in 1:2:n-1
                push!(gates,B.Rotation(B.word([:Z,:Z],[q,q+1]),0.37))
            end
            if family==:weak
                for j in 1:coupling
                    q=2mod(j-1,max(1,div(n,2)-1))+2
                    q<n && push!(gates,B.Rotation(B.word([:Z,:Z],[q,q+1]),0.19))
                end
            elseif family==:cross
                for _ in 1:coupling*max(1,div(n,2))
                    a=rand(rng,1:div(n,2)); b=rand(rng,div(n,2)+1:n)
                    push!(gates,B.Clifford(:CNOT,a,b))
                end
            end
        elseif family==:ising
            for q in 1:n-1
                push!(gates,B.Rotation(B.word([:Z,:Z],[q,q+1]),0.31))
            end
            periodic && push!(gates,B.Rotation(B.word([:Z,:Z],[n,1]),0.31))
            for q in 1:n
                push!(gates,B.Rotation(B.word([:X],[q]),0.47))
                push!(gates,B.Rotation(B.word([:Z],[q]),0.23))
            end
        elseif family==:clifford
            for q in 1:n
                push!(gates,B.Clifford(isodd(q+layer) ? :H : :S,q))
            end
            for q in 1:2:n-1; push!(gates,B.Clifford(:CNOT,q,q+1)); end
        end
    end
    B.Circuit(n,gates)
end
function observable(n,kind)
    if kind==:single
        return [B.word([:Z],[cld(n,2)])=>1.0]
    elseif kind==:zz
        return [B.word([:Z,:Z],[cld(n,2),cld(n,2)+1])=>1.0]
    elseif kind==:local_sum
        return [B.word([:Z],[q])=>1/n for q in 1:n]
    else
        return [B.word([:X,:X],[1,n])=>1.0]
    end
end
function runone(method,c,terms; measure_memory=false)
    t0=time_ns(); init=0; prop=0; agg=0; peak=0; op=nothing
    if method=="S0"
        d=Dict{PW,Float64}()
        for (p,v) in terms; d[p]=get(d,p,0.0)+v; end
        init=time_ns()-t0; peak=measure_memory ? Base.summarysize(d) : 0
        t=time_ns()
        for g in Iterators.reverse(c.gates)
            next=Dict{PW,Float64}()
            for (p,v) in d, (q,w) in B.expand(g,p,v)
                w==0 && continue
                next[q]=get(next,q,0.0)+w
                length(next)>lim.max_contributions && throw(ResourceLimit("sequential contribution limit exceeded"))
            end
            d=Dict(p=>v for (p,v) in next if v!=0.0)
            measure_memory && (peak=max(peak,Base.summarysize(d)))
        end
        prop=time_ns()-t
        return (init,prop,agg,time_ns()-t0,peak,0,0,0,0,0,expectation_zero(d),d)
    end
    leaf=method=="S1" ? cld(c.nqubits,2) : 2
    op=Operator(c.nqubits,leaf,terms;limits=lim)
    init=time_ns()-t0; peak=measure_memory ? Base.summarysize(op) : 0
    t=time_ns()
    for g in Iterators.reverse(c.gates)
        SpatialPauli.propagate!(op,g)
        measure_memory && (peak=max(peak,Base.summarysize(op)))
    end
    prop=time_ns()-t; agg=op.aggregation_ns
    s=stats(op)
    (init,prop,agg,time_ns()-t0,peak,s.peak_components,s.peak_terms,s.peak_contributions,s.cross_gates,s.before_aggregation,expectation_zero(op),op)
end
function csvline(xs)
    join((replace(string(x),'"'=>"''",','=>';') for x in xs),",")
end
function main()
    root=joinpath(@__DIR__,"..","results","raw"); mkpath(root)
    path=joinpath(root,"benchmark_"*Dates.format(now(),"yyyymmdd_HHMMSS")*".csv")
    header="family,n,depth,observable,coupling,periodic,method,status,median_total_s,min_total_s,max_total_s,median_init_s,median_propagate_s,median_aggregate_s,peak_live_bytes,allocated_bytes,peak_components,peak_local_terms,peak_contributions,cross_gates,before_aggregation,expectation,error,Julia,OS,threads,CPU,seed,message"
    open(path,"w") do io; println(io,header); end
    ns=parse.(Int,split(get(ENV,"SPATIAL_NS","4,8,12,16"),","))
    depths=parse.(Int,split(get(ENV,"SPATIAL_DEPTHS","1,2"),","))
    sweep=get(ENV,"SPATIAL_COUPLING_SWEEP","")=="1"
    families=sweep ? [:weak] : Symbol.(split(get(ENV,"SPATIAL_FAMILIES","pairs,ising,weak,cross,clifford"),","))
    kinds=sweep ? (:local_sum,) : Tuple(Symbol.(split(get(ENV,"SPATIAL_OBSERVABLES","single,zz,local_sum,multi_block"),",")))
    println("Writing ",path)
    for n in ns, depth in depths, family in families, kind in kinds, coupling_override in (sweep ? (0:3) : (-1:-1))
        seed=1729+n*31+depth*7+Int(findfirst(==(family),families))
        coupling=sweep ? coupling_override : (family==:weak || family==:cross ? depth : 0)
        c=circuit(n,family,depth,seed;coupling)
        terms=observable(n,kind)
        ref=nothing
        try
            ref=sequential(c,terms;limits=lim)
        catch e
            # Method rows still record independent statuses below.
        end
        for method in ("S0","S1","S2")
            status="ok"; msg=""; samples=Any[]; bytes=0; memory_result=nothing
            try
                runone(method,c,terms) # warmup, excluded
                for _ in 1:3
                    GC.gc()
                    b=@allocated result=runone(method,c,terms)
                    bytes=max(bytes,b); push!(samples,result)
                end
                memory_result=runone(method,c,terms;measure_memory=true)
            catch e
                status=e isa ResourceLimit ? "resource_limited" : "error"
                msg=sprint(showerror,e)
            end
            vals=Any[family,n,depth,kind,coupling,false,method,status]
            if status=="ok"
                totals=[x[4]/1e9 for x in samples]; k=sortperm(totals)[2]; r=samples[k]
                err=ref===nothing ? "" : abs(r[11]-expectation_zero(ref))
                append!(vals,[median(totals),minimum(totals),maximum(totals),r[1]/1e9,r[2]/1e9,r[3]/1e9,memory_result[5],bytes,r[6],r[7],r[8],r[9],r[10],r[11],err])
            else
                append!(vals,fill("",15))
            end
            append!(vals,[VERSION,Sys.KERNEL,Threads.nthreads(),Sys.CPU_NAME,seed,msg])
            open(path,"a") do io; println(io,csvline(vals)); end
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
