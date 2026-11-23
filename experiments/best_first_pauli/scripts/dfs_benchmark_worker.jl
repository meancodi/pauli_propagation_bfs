using JSON3, TOML, Dates, LinearAlgebra
include(joinpath(@__DIR__,"..","src","BestFirstPauli.jl"))
using .BestFirstPauli
const B=BestFirstPauli
const PP=B.PP
LinearAlgebra.BLAS.set_num_threads(1)
Threads.nthreads()==1 || error("DFS comparison requires --threads=1")
const ROOT=normpath(joinpath(@__DIR__,".."))
const OUTROOT=joinpath(ROOT,"results_dfs_vs_bestfirst")
cfg=TOML.parsefile(joinpath(ROOT,"dfs_benchmark_config.toml"))
n=parse(Int,ARGS[1]); skip=parse(Int,ARGS[2]); outfile,progressfile=ARGS[3:4]
cases=[(family,layer,seed,obs) for family in cfg["families"] for layer in cfg["layers"]
       for seed in cfg["seeds"] for obs in cfg["observables"] if seed==11 || obs=="Z_center"]
function progress(i,phase)
    open(progressfile,"a") do io
        println(io,JSON3.write((;case_index=i,phase,timestamp=string(now(UTC)))))
    end
end
function appendrow(row)
    open(outfile,"a") do io
        println(io,JSON3.write(row))
    end
end
function search(algorithm,c,p;kwargs...)
    algorithm=="bestfirst" ? bestfirst(c,p;kwargs...) : depthfirst(c,p;kwargs...)
end
function baserow(i,c,family,layer,seed,obs,algorithm,experiment,budget,rep,reference,reference_ms,total_rotations)
    Dict{String,Any}(
      "run_id"=>"$(n)_$(family)_$(layer)_$(seed)_$(obs)_$(algorithm)_$(experiment)_$(budget)_$(rep)",
      "case_index"=>i,"algorithm"=>algorithm,"nqubits"=>n,"family"=>family,"layers"=>layer,
      "seed"=>seed,"observable"=>obs,"experiment"=>experiment,"requested_budget"=>budget,
      "repeat"=>rep,"reference"=>reference,"reference_backend"=>"VectorPauliSum exact",
      "reference_elapsed_ms"=>reference_ms,"total_gates"=>length(c.gates),
      "total_rotations"=>total_rotations,"threads"=>Threads.nthreads(),"julia_version"=>string(VERSION),
      "timestamp"=>string(now(UTC))*"Z")
end
function record!(row,timed,r)
    row["status"]=string(r.status); row["reason"]=string(r.reason); row["expectation"]=r.expectation
    row["absolute_error"]=ismissing(row["reference"]) ? missing : abs(r.expectation-row["reference"])
    row["relative_error"]=ismissing(row["reference"]) || abs(row["reference"])<=1e-12 ? missing : row["absolute_error"]/abs(row["reference"])
    row["frontier_bound"]=r.frontier_bound; row["requested_accuracy_achieved"]=
      row["experiment"]=="accuracy_tolerance" ? r.status==:certified || r.status==:exact : missing
    row["elapsed_ms"]=timed.time*1000; row["search_elapsed_ms"]=r.elapsed_seconds*1000
    row["allocated_bytes"]=timed.bytes; row["gc_ms"]=timed.gctime*1000
    row["peak_rss_mib"]=Sys.maxrss()/2.0^20
    row["rss_source"]="Sys.maxrss process-lifetime high-water; includes JIT/warmup/prior cases"
    for key in (:expanded_nodes,:generated_children,:completed_leaves,:peak_frontier,:frontier_nodes,
                :queue_pushes,:queue_pops,:stack_pushes,:stack_pops,:mean_frontier,
                :approx_node_bytes,:max_tree_depth,:bookkeeping_seconds,:kernel_seconds)
        row[string(key)]=getproperty(r,key)
    end
    row["branch_events_lower_bound"]=max(0,r.generated_children-r.expanded_nodes)
    row["observed_children_per_expansion"]=r.expanded_nodes==0 ? missing : r.generated_children/r.expanded_nodes
    appendrow(row)
end
for (i,(family,layer,seed,obs)) in enumerate(cases)
    i<=skip && continue
    progress(i,"reference")
    c=make_circuit(n,Symbol(family),layer,seed); p=observable(n,Symbol(obs))
    rotations=count(g->g isa Rotation,c.gates)
    reference=missing; reference_ms=missing
    try
        gates,theta,ob=B.upstream_inputs(c,p); init=PP.VectorPauliSum(ob)
        PP.propagate(gates,init,theta;min_abs_coeff=0.0,max_weight=Inf,thread=false) # warmup
        t=@timed PP.propagate(gates,init,theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
        reference=PP.overlapwithzero(t.value); reference_ms=t.time*1000
    catch err
        appendrow(Dict("run_id"=>"reference_failure_$(n)_$(i)","case_index"=>i,"algorithm"=>"upstream",
          "nqubits"=>n,"family"=>family,"layers"=>layer,"seed"=>seed,"observable"=>obs,
          "experiment"=>"reference","status"=>"failed","reason"=>sprint(showerror,err),
          "timestamp"=>string(now(UTC))*"Z"))
    end
    # Compile both paths and all keyword variants outside measurements.
    search("bestfirst",c,p;max_expansions=10,max_frontier_nodes=cfg["max_frontier_nodes"],max_runtime_seconds=Inf)
    search("dfs",c,p;max_expansions=10,max_frontier_nodes=cfg["max_frontier_nodes"],max_runtime_seconds=Inf)
    experiments=Tuple{String,Any}[("expansion_budget",b) for b in cfg["expansion_budgets"]]
    append!(experiments,[("runtime_budget",b) for b in cfg["runtime_budgets_seconds"]])
    append!(experiments,[("accuracy_tolerance",b) for b in cfg["accuracy_tolerances"]])
    for (experiment,budget) in experiments
        for rep in 1:cfg["repeats"]
            # Alternate order by case/repeat/budget, reducing systematic thermal/order bias.
            algorithms=isodd(i+rep+findfirst(==(budget),last.(experiments))) ? ("dfs","bestfirst") : ("bestfirst","dfs")
            for algorithm in algorithms
                GC.gc(); progress(i,"$(experiment)_$(budget)_$(rep)_$(algorithm)")
                kwargs = experiment=="expansion_budget" ?
                  (;max_expansions=Int(budget),max_frontier_nodes=cfg["max_frontier_nodes"],max_runtime_seconds=Inf) :
                  experiment=="runtime_budget" ?
                  (;max_expansions=typemax(Int),max_frontier_nodes=cfg["max_frontier_nodes"],max_runtime_seconds=Float64(budget)) :
                  (;max_expansions=cfg["accuracy_max_expansions"],max_frontier_nodes=cfg["max_frontier_nodes"],
                    max_runtime_seconds=cfg["accuracy_max_runtime_seconds"],error_tolerance=Float64(budget))
                timed=@timed search(algorithm,c,p;kwargs...)
                row=baserow(i,c,family,layer,seed,obs,algorithm,experiment,budget,rep,reference,reference_ms,rotations)
                record!(row,timed,timed.value)
            end
        end
    end
    # Separate instrumented run: useful for attribution, excluded from primary timings.
    for algorithm in ("dfs","bestfirst")
        GC.gc(); progress(i,"overhead_$(algorithm)")
        timed=@timed search(algorithm,c,p;max_expansions=cfg["overhead_expansions"],
          max_frontier_nodes=cfg["max_frontier_nodes"],max_runtime_seconds=Inf,measure_bookkeeping=true)
        row=baserow(i,c,family,layer,seed,obs,algorithm,"overhead_profile",cfg["overhead_expansions"],1,reference,reference_ms,rotations)
        record!(row,timed,timed.value)
    end
    println("n=$n DFS/BFS case $i/$(length(cases)) $family L=$layer seed=$seed $obs done"); flush(stdout)
end
progress(length(cases)+1,"done")
