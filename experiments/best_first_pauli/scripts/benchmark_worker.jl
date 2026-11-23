using JSON3, TOML, Dates, LinearAlgebra
include(joinpath(@__DIR__,"..","src","BestFirstPauli.jl"))
using .BestFirstPauli
const B = BestFirstPauli
const PP = B.PP
LinearAlgebra.BLAS.set_num_threads(1)
Threads.nthreads() == 1 || error("Primary suite requires --threads=1")
const ROOT = normpath(joinpath(@__DIR__,".."))
config = TOML.parsefile(joinpath(ROOT,"benchmark_config.toml"))
n = parse(Int,ARGS[1]); backend = ARGS[2]; skip = parse(Int,ARGS[3])
outfile,progressfile = ARGS[4:5]
manifest = TOML.parsefile(joinpath(ROOT,"Manifest.toml"))
pkg = only(manifest["deps"]["PauliPropagation"])
cases = [(family,layer,seed,obs,0.0,"primary") for family in config["families"] for layer in config["layers"] for seed in config["seeds"] for obs in config["observables"]]
if backend == "bestfirst"
    append!(cases,[(family,4,11,"ZZ_center",budget,"budget") for family in config["families"] for budget in config["budget_seconds"]])
end
function progress(i,phase)
    # Append-only events avoid Windows rename/delete sharing races with the monitor.
    open(progressfile,"a") do io
        println(io,JSON3.write((;case_index=i,phase,timestamp=string(now(UTC)))))
    end
end
function measure_bfs(c,p,config,cap)
    bestfirst(c,p; max_expansions=config["max_expansions"],max_frontier_nodes=config["max_frontier_nodes"],max_runtime_seconds=cap)
end
for (i,(family,layer,seed,obs,budget,suite)) in enumerate(cases)
    i <= skip && continue
    progress(i,"preparation")
    c = make_circuit(n,Symbol(family),layer,seed); p = observable(n,Symbol(obs))
    cap = suite == "budget" ? budget : config["max_runtime_seconds"]
    gates,theta,ob = B.upstream_inputs(c,p)
    # All representation conversion and angle preparation is outside timing.
    init = backend == "VectorPauliSum" ? PP.VectorPauliSum(ob) : PP.PauliSum(ob)
    runone = backend == "bestfirst" ? ()->measure_bfs(c,p,config,cap) : ()->begin
        out=PP.propagate(gates,init,theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
        (;expectation=PP.overlapwithzero(out),unique_terms=length(out))
    end
    progress(i,"warmup")
    runone()
    peakterms = backend == "bestfirst" ? missing : PP.@peakpaulis runone()
    dense = n == 5 ? dense_expectation(c,p) : missing
    for rep in 1:config["repeats"]
        GC.gc()
        progress(i,"repeat_$rep")
        timed = @timed runone()
        r = timed.value
        row = Dict{String,Any}(
            "run_id"=>"$(n)_$(family)_$(layer)_$(seed)_$(obs)_$(backend)_$(suite)_$(budget)_$(rep)",
            "algorithm"=>backend == "bestfirst" ? "best_first" : "upstream", "backend"=>backend,"case_index"=>i,
            "nqubits"=>n,"family"=>family,"layers"=>layer,"seed"=>seed,"observable"=>obs,
            "theta_spec"=>JSON3.write(theta),"truncation_config"=>"min_abs_coeff=0;max_weight=Inf",
            "priority_mode"=>backend == "bestfirst" ? "abs_coeff" : missing,"suite"=>suite,"repeat"=>rep,
            "status"=>backend == "bestfirst" ? string(r.status) : "exact",
            "reason"=>backend == "bestfirst" ? string(r.reason) : "completed",
            "expectation"=>r.expectation,
            "dense_reference"=>dense,"elapsed_ms"=>timed.time*1000,"allocated_bytes"=>timed.bytes,
            "gc_ms"=>timed.gctime*1000,"peak_rss_mib"=>Sys.maxrss()/2.0^20,
            "rss_source"=>"Sys.maxrss: process lifetime high-water, includes compilation/warmup/prior cases",
            "peak_terms"=>peakterms,"unique_terms"=>backend == "bestfirst" ? missing : r.unique_terms,
            "max_expansions"=>backend == "bestfirst" ? config["max_expansions"] : missing,
            "max_frontier_nodes"=>backend == "bestfirst" ? config["max_frontier_nodes"] : missing,
            "runtime_cap_seconds"=>backend == "bestfirst" ? cap : missing,
            "julia_version"=>string(VERSION),"package_revision"=>pkg["git-tree-sha1"],
            "package_version"=>pkg["version"],"threads"=>Threads.nthreads(),"timestamp"=>string(now(UTC))*"Z")
        for key in (:frontier_bound,:frontier_abs_coeff_sum,:expanded_nodes,:generated_children,:completed_leaves,
                    :peak_frontier,:queue_pushes,:queue_pops,:max_tree_depth,:frontier_nodes,:approx_node_bytes,:mean_frontier)
            row[string(key)] = backend == "bestfirst" ? getproperty(r,key) : missing
        end
        open(outfile,"a") do io
            println(io,JSON3.write(row))
        end
    end
    # Exact parameter/circuit specification independent of RNG implementation changes.
    circuitfile = joinpath(ROOT,"results","raw","circuits","$(n)_$(family)_$(layer)_$(seed).json")
    if !isfile(circuitfile)
        spec = [g isa Rotation ? (;kind="rotation",symbols=symbols(g.generator,n),theta=g.theta) : (;kind=string(g.kind),sites=[g.a,g.b]) for g in c.gates]
        write(circuitfile,JSON3.write((;nqubits=n,family,layer,seed,gates=spec)))
    end
    println("n=$n $backend case $i/$(length(cases)) $family L=$layer $obs $suite done"); flush(stdout)
end
progress(length(cases)+1,"done")
