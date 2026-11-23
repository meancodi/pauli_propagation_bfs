using JSON3, LinearAlgebra, Statistics, Random, Dates, SHA
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","FavorableCircuits.jl")); using .FavorableCircuits
include(joinpath(@__DIR__,"..","..","hierarchical_pauli","src","HierarchicalPauli.jl")); using .HierarchicalPauli
include(joinpath(@__DIR__,"..","..","hierarchical_pauli","src","HierarchicalPauliV2.jl")); using .HierarchicalPauliV2
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); outdir=joinpath(root,"results","favorable_v1"); mkpath(outdir); out=joinpath(outdir,"campaign.jsonl"); progress=joinpath(outdir,"campaign_progress.log")
lock=joinpath(outdir,"campaign.lock"); force_lock=any(==("--force-lock"),ARGS)
isdir(lock) && (force_lock ? rm(lock;recursive=true,force=true) : error("campaign lock exists"))
mkdir(lock); atexit(()->isdir(lock)&&rm(lock;recursive=true,force=true))
smoke=any(==("--smoke"),ARGS); repeats=smoke ? 3 : 5; done=Set{String}(); isfile(out)&&for l in eachline(out); try push!(done,String(JSON3.read(l).config_id)) catch end end
function logline(x); open(progress,"a") do io; println(io,Dates.format(now(),"yyyy-mm-dd HH:MM:SS")," ",x); end end
function emit(x); open(out,"a") do io; println(io,JSON3.write(x)); end end
function applysum(m,d)
    o=Dict{PauliWord,Float64}(); for (p,c) in d; for (q,v) in Happly(m,p,c); o[q]=get(o,q,0.0)+v end end; HierarchicalPauli.merge_terms(o)
end
Happly(m,p,c)=HierarchicalPauli.apply_map(m,p,c)
function cached_local(c,p,b)
    d=Dict(p=>1.0); cache=Dict{UInt64,Any}(); built=0.0; chunks=0
    for last in length(c.gates):-b:1
        first=max(1,last-b+1); gs=c.gates[first:last]; key=UInt64(hash(gs)); m=get(cache,key,nothing)
        if m===nothing; z=@timed HierarchicalPauli.build_exact_block_map(gs); m=z.value; cache[key]=m; built+=z.time end
        d=applysum(m,d); chunks+=1
    end
    (;expectation=sum((v for (q,v) in d if q.x==0);init=0.0),terms=length(d),built_seconds=built,cache_size=length(cache),cache_hits=chunks-length(cache))
end
function repeated_cached(block,p,reps;powers=false)
    z=@timed HierarchicalPauli.build_exact_block_map(block); base=z.value; powerset=Any[base]; power_build=z.time; maxpow=1
    if powers
        while 2maxpow<=reps
            q=HierarchicalPauli.compose_maps(powerset[end],powerset[end]);
            sum(length(v) for v in values(q.columns))>5000 && break
            push!(powerset,q); maxpow*=2; power_build+=q.compose_seconds
        end
    end
    d=Dict(p=>1.0); applied=0; i=length(powerset)
    while applied<reps
        take=min(2^(i-1),reps-applied); d=applysum(powerset[i],d); applied+=take
        i==1 && continue; take<2^(i-1) && (i-=1)
    end
    (;expectation=sum((v for (q,v) in d if q.x==0);init=0.0),terms=length(d),built_seconds=z.time,power_seconds=power_build,materialized_powers=length(powerset),applied)
end
function obslist(c)
    o=FavorableCircuits.observables(c.nqubits); [o.X,o.Z,o.ZZ]
end
function make_case(family, cfg)
    family==:echo && return echo_circuit(cfg.n,cfg.depth,cfg.seed;perturb=cfg.perturb,theta=cfg.theta)
    family==:repeated && return repeated_block_circuit(cfg.k,cfg.ngates,cfg.reps,cfg.seed;kind=cfg.kind,perturbed=cfg.vary)
    family==:clifford_sandwich && return clifford_sandwich_circuit(cfg.n,cfg.depth,cfg.nc,cfg.seed;theta=cfg.theta)
    family==:trotter && return trotter_circuit(cfg.n,cfg.reps,cfg.dt,cfg.perturb;vary=cfg.vary,seed=cfg.seed)
    periodic_circuit(cfg.n,cfg.reps,cfg.perturb;order=cfg.order,seed=cfg.seed)
end
function run_method(method,c,p,cfg)
    th=1e-3
    method==:B0 && return upstream(c,p,:VectorPauliSum)
    method==:B2 && return propagate(c,p;block_size=1,threshold=th)
    method==:B3 && return propagate(c,p;block_size=cfg.block_size,threshold=th)
    method==:B4 && return cached_local(c,p,cfg.cache_block)
    method==:B5 && return v2_propagate(c.gates,p;strategy=:balanced,config=make_config(truncation_mode=:all_internal,coefficient_threshold=th,max_support=4,max_columns=256))
    method==:B6 && return v2_propagate(c.gates,p;strategy=:support_constrained,config=make_config(truncation_mode=:all_internal,coefficient_threshold=th,max_support=3,max_columns=64))
    method==:B8 && return v2_propagate(c.gates,p;strategy=:balanced,config=make_config(max_support=0,max_columns=1))
    method==:B7 && return repeated_cached(cfg.block,c,cfg.reps;powers=true) # corrected below for repeated cases
end
function runone(method,c,p,cfg)
    if method==:B7
        hasproperty(cfg,:block) || return (;status=:not_applicable,expectation=missing)
        return repeated_cached(cfg.block,p,cfg.reps;powers=true)
    end
    run_method(method,c,p,cfg)
end
methods=(:B0,:B2,:B3,:B4,:B5,:B6,:B7,:B8)
logline("START smoke=$smoke repeats=$repeats")
for family in (:echo,:repeated,:clifford_sandwich,:trotter,:periodic)
    cfgs=if smoke
        family==:echo ? [(;n=3,depth=4,seed=11,perturb=:noncommuting,theta=.1,block_size=4,cache_block=4)] :
        family==:repeated ? [(;k=3,ngates=4,reps=4,seed=11,kind=:mixed,vary=false,block=supported_local_block(3,4,11;kind=:mixed),block_size=4,cache_block=4)] :
        family==:clifford_sandwich ? [(;n=4,depth=4,nc=2,seed=11,theta=.2,block_size=4,cache_block=4)] :
        family==:trotter ? [(;n=4,reps=4,dt=.05,perturb=.01,vary=false,seed=11,block_size=4,cache_block=4)] :
        [(;n=3,reps=8,perturb=.01,order=4,seed=11,block_size=4,cache_block=4)]
    else
        family==:echo ? [(;n=n,depth=d,seed=s,perturb=pt,theta=th,block_size=4,cache_block=4) for n in (3,4),d in (4,16),s in (11,29),pt in (:identity,:noncommuting),th in (.05,.2)][1:16] :
        family==:repeated ? [(;k=k,ngates=g,reps=r,seed=11,kind=kind,vary=false,block=supported_local_block(k,g,11;kind=kind),block_size=g,cache_block=g) for k in (2,3,4),g in (4,8),r in (2,8,32),kind in (:clifford,:mixed)] :
        family==:clifford_sandwich ? [(;n=n,depth=d,nc=nc,seed=11,theta=th,block_size=4,cache_block=4) for n in (4,8),d in (4,16),nc in (1,4),th in (.1,.3)] :
        family==:trotter ? [(;n=n,reps=r,dt=dt,perturb=pt,vary=false,seed=11,block_size=4,cache_block=4) for n in (4,8),r in (2,8,16),dt in (.01,.05),pt in (0.0,.05)] :
        [(;n=n,reps=r,perturb=pt,order=4,seed=11,block_size=4,cache_block=4) for n in (2,3,4),r in (4,16,32),pt in (0.0,.01,.05)]
    end
    for cfg in cfgs
        c=make_case(family,cfg); for p in obslist(c), method in methods
            stable=bytes2hex(sha1(codeunits(string(family,"|",cfg,"|",p.x,"|",p.z,"|",method))))
            id=string(family,"_",stable)
            id in done && continue
            ref=try upstream(c,p,:VectorPauliSum).expectation catch; missing end
            try runone(method,c,p,cfg); runone(method,c,p,cfg) catch end
            for rep in 1:repeats
                z=@timed begin
                    try; runone(method,c,p,cfg) catch ex; (;status=:failed,reason=sprint(showerror,ex),expectation=missing) end
                end
                r=z.value; ex=hasproperty(r,:expectation) ? r.expectation : missing
                emit((;config_id=id*"_"*string(rep),family=String(family),method=String(method),rep,qubits=c.nqubits,gates=length(c.gates),reference=ref,expectation=ex,observed_error=ismissing(ex)||ismissing(ref) ? missing : abs(ex-ref),seconds=z.time,allocated_bytes=z.bytes,status=hasproperty(r,:status) ? String(r.status) : "complete",terms=hasproperty(r,:terms) ? length(r.terms) : missing,observables_x=string(p.x),observables_z=string(p.z),timestamp=string(now(UTC))))
            end
            push!(done,id); logline("done $id")
        end
    end
end
logline("FINISH")
