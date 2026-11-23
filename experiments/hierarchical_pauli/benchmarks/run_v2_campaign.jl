using JSON3, Random, Dates, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","HierarchicalPauli.jl")); using .HierarchicalPauli
include(joinpath(@__DIR__,"..","src","HierarchicalPauliV2.jl")); using .HierarchicalPauliV2
BLAS.set_num_threads(1)

root=normpath(joinpath(@__DIR__,"..")); outdir=joinpath(root,"results","hierarchical_v2"); mkpath(outdir)
lockdir=joinpath(outdir,"campaign.lock")
isdir(lockdir) && error("campaign lock exists at $lockdir; do not launch a second runner")
mkdir(lockdir); open(joinpath(lockdir,"pid.txt"),"w") do io; print(io,getpid()) end
atexit(() -> isdir(lockdir) && rm(lockdir;recursive=true,force=true))
outfile=joinpath(outdir,"benchmark.jsonl"); progress=joinpath(outdir,"overnight_progress.log"); manifest=joinpath(outdir,"run_manifest.json")
tier=any(==("--smoke"),ARGS) ? :smoke : (any(==("--stress"),ARGS) ? :stress : :moderate)
ns=tier==:smoke ? [5,8] : tier==:moderate ? [5,8,12,16] : [8,12,16,24]
layers=tier==:smoke ? [1,2] : tier==:moderate ? [1,2,4] : [2,4,8]
families=[:ising,:random]; observables=[:X_center,:Y_center,:Z_center,:ZZ_center]; seeds=[11,29]
repeats=tier==:smoke ? 2 : 3; started=time(); wall_budget=tier==:stress ? 8*3600.0 : 2*3600.0
completed=Set{String}(); isfile(outfile) && for line in eachline(outfile); try push!(completed,String(JSON3.read(line).config_id)) catch end end
function cancellation_circuit(n,layers,seed)
    rng=MersenneTwister(seed); gs=SupportedGate[]
    for _ in 1:layers, q in 1:n-1
        p=word([:Z,:Z],[q,q+1]); push!(gs,Rotation(p,.42)); push!(gs,Rotation(p,-.42+.015rand(rng)))
    end
    Circuit(n,gs)
end
function makeobs(n,k)
    k==:X_center && return word([:X],[cld(n,2)]); k==:Y_center && return word([:Y],[cld(n,2)])
    k==:Z_center && return word([:Z],[cld(n,2)]); k==:ZZ_center && return word([:Z,:Z],[cld(n,2),cld(n,2)+1])
end
function emit(row)
    open(outfile,"a") do io; println(io,JSON3.write(row)); end
end
function logline(s)
    open(progress,"a") do io; println(io,Dates.format(now(),"yyyy-mm-dd HH:MM:SS")," ",s); end
end
function timed_run(f)
    GC.gc(); z=@timed f(); (z.value,z.time,z.bytes,z.gctime)
end
function run_method(name,c,p,threshold)
    if name==:B0; r=upstream(c,p,:VectorPauliSum); return r, :complete
    elseif name==:B1; return propagate(c,p;block_size=1,threshold=threshold), :complete
    elseif name==:B2; return propagate(c,p;block_size=4,threshold=threshold), :complete
    elseif name==:B3; return v2_propagate(c.gates,p;strategy=:balanced,config=make_config(truncation_mode=:none,max_support=4,max_columns=256)), :v2
    elseif name==:B4; return v2_propagate(c.gates,p;strategy=:balanced,config=make_config(truncation_mode=:all_internal,coefficient_threshold=threshold,max_support=4,max_columns=256)), :v2
    elseif name==:B5; return v2_propagate(c.gates,p;strategy=:left_deep,config=make_config(truncation_mode=:all_internal,coefficient_threshold=threshold,max_support=4,max_columns=256)), :v2
    elseif name==:B6; return v2_propagate(c.gates,p;strategy=:support_constrained,config=make_config(truncation_mode=:all_internal,coefficient_threshold=threshold,max_support=3,max_columns=64)), :v2
    elseif name==:B7; return v2_propagate(c.gates,p;strategy=:balanced,config=make_config(truncation_mode=:none,max_support=1,max_columns=1)), :v2
    end
end
methods=(:B0,:B1,:B2,:B3,:B4,:B5,:B6,:B7)
logline("START tier=$tier repeats=$repeats")
for n in ns, l in layers, family in families, ok in observables, seed in seeds
    time()-started>wall_budget && break
    c=make_circuit(n,family,l,seed); p=try makeobs(n,ok) catch; continue end
    # Use the same exact upstream reference for every method. Dense is also
    # recorded for the small cases but is never substituted for a missing ref.
    reference=try upstream(c,p,:VectorPauliSum).expectation catch; missing end
    for threshold in (1e-3,1e-2), method in methods, rep in 1:repeats
        id="$(tier)_$(n)_$(l)_$(family)_$(ok)_$(seed)_$(method)_$(threshold)_$(rep)"
        id in completed && continue
        t0=time(); status="complete"; reason="finished"; value=missing; elapsed=missing; bytes=missing; gc=missing
        try
            f=()->run_method(method,c,p,threshold); f(); value,elapsed,bytes,gc=timed_run(f)
            r=value[1]; status=hasproperty(r,:status) ? String(r.status) : "complete"
            reason=hasproperty(r,:reason) ? String(r.reason) : "finished"
            obs=hasproperty(r,:expectation) ? r.expectation : missing
            mm=hasproperty(r,:metrics) ? r.metrics : nothing
            emit((;config_id=id,tier=String(tier),nqubits=n,layers=l,gates=length(c.gates),family=String(family),observable=String(ok),seed,method=String(method),threshold,repeat=rep,reference,expectation=obs,observed_error=ismissing(reference)||ismissing(obs) ? missing : abs(obs-reference),elapsed_seconds=elapsed,allocated_bytes=bytes,gc_seconds=gc,status,reason,peak_terms=hasproperty(r,:terms) ? length(r.terms) : missing,certified_bound=hasproperty(r,:certified_bound) ? r.certified_bound : missing,discarded_mass=mm===nothing ? missing : mm.discarded_mass,internal_truncations=mm===nothing ? missing : mm.internal_truncations,materialized_nodes=mm===nothing ? missing : mm.materialized_nodes,sequential_nodes=mm===nothing ? missing : mm.sequential_nodes,resource_events=mm===nothing ? missing : mm.resource_events,cache_hits=mm===nothing ? missing : mm.cache_hits,cache_misses=mm===nothing ? missing : mm.cache_misses,tree_seconds=mm===nothing ? missing : mm.tree_seconds,composition_seconds=mm===nothing ? missing : mm.composition_seconds,package_version="PauliPropagation 0.9.2",julia_version=string(VERSION),threads=Threads.nthreads(),timestamp=string(now(UTC))))
            push!(completed,id); logline("done $id runtime=$(elapsed)")
        catch ex
            emit((;config_id=id,tier=String(tier),nqubits=n,layers=l,gates=length(c.gates),family=String(family),observable=String(ok),seed,method=String(method),threshold,repeat=rep,reference,status="failed",reason=sprint(showerror,ex),elapsed_seconds=time()-t0,package_version="PauliPropagation 0.9.2",julia_version=string(VERSION)))
            logline("failed $id $(sprint(showerror,ex))")
        end
    end
end
open(manifest,"w") do io; JSON3.write(io,(;tier=String(tier),started=string(now(UTC)),julia_version=string(VERSION),threads=Threads.nthreads(),wall_budget_seconds=wall_budget,resumable=true,hard_memory_limit=false,warning="Windows runner uses admission estimates and per-case Julia safeguards; no hard RSS kill is claimed.")); end
logline("FINISH tier=$tier")
