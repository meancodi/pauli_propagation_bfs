using JSON3,TOML,Dates,SHA,LinearAlgebra
include(joinpath(@__DIR__,"..","src","telemetry.jl"))
const ROOT=normpath(joinpath(@__DIR__,"..")); const OUTROOT=joinpath(ROOT,"results_dfs_vs_bestfirst")
Threads.nthreads()==1 || error("Run with --threads=1"); LinearAlgebra.BLAS.set_num_threads(1)
cfg=TOML.parsefile(joinpath(ROOT,"dfs_benchmark_config.toml"))
selected=isempty(ARGS) ? cfg["nqubits"] : parse.(Int,ARGS)
all(in(cfg["nqubits"]),selected) || error("Supported sizes: $(cfg["nqubits"])")
for p in ("raw","work","figures"); mkpath(joinpath(OUTROOT,p)); end
metadata=Dict("config"=>cfg,"julia_version"=>string(VERSION),"os"=>string(Sys.KERNEL),
 "os_version"=>Sys.iswindows() ? strip(read(`cmd /c ver`,String)) : strip(read(`uname -a`,String)),
 "machine"=>Sys.MACHINE,"cpu"=>Sys.cpu_info()[1].model,"logical_cpu_count"=>Sys.CPU_THREADS,
 "total_memory_bytes"=>Sys.total_memory(),"threads"=>1,"blas_threads"=>1,"timestamp"=>string(now(UTC))*"Z",
 "gpu"=>try strip(read(`nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader`,String)) catch; "not detected; unused" end,
 "source_sha256"=>Dict(relpath(joinpath(d,f),ROOT)=>bytes2hex(sha256(read(joinpath(d,f)))) for (d,_,fs) in walkdir(joinpath(ROOT,"src")) for f in fs),
 "manifest_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml")))),"rss_poll_interval_seconds"=>0.05)
metafile=joinpath(OUTROOT,"raw","metadata.json")
if isfile(metafile)
    old=JSON3.read(read(metafile,String),Dict{String,Any})
    for k in ("config","source_sha256","manifest_sha256","julia_version")
        old[k]==metadata[k] || error("Existing extension results use different $k; archive results_dfs_vs_bestfirst first")
    end
else
    write(metafile,JSON3.write(metadata))
end
for n in selected
    outfile=joinpath(OUTROOT,"raw","$(n)_search.jsonl"); progressfile=joinpath(OUTROOT,"work","progress_$(n).jsonl")
    skip=0
    if isfile(outfile)
        prior=[JSON3.read(s) for s in eachline(outfile)]
        completed=[Int(r.case_index) for r in prior if get(r,:experiment,"")=="overhead_profile" && get(r,:algorithm,"")=="bestfirst"]
        skip=maximum(completed;init=0)
    end
    while true
        write(progressfile,JSON3.write((;case_index=skip+1,phase="startup"))*"\n")
        logfile=joinpath(OUTROOT,"work","$(n)_after_$(skip).log"); logio=open(logfile,"w")
        cmd=`$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$ROOT $(joinpath(@__DIR__,"dfs_benchmark_worker.jl")) $n $skip $outfile $progressfile`
        proc=run(pipeline(cmd;stdout=logio,stderr=logio);wait=false); pid=getpid(proc)
        current=skip+1; began=time(); peak=0; stopreason=""
        while process_running(proc)
            sleep(.05)
            state=try JSON3.read(last(readlines(progressfile))) catch; nothing end
            if state!==nothing && state.case_index!=current; current=state.case_index; began=time(); end
            rss=process_rss(pid)
            if !ismissing(rss)
                peak=max(peak,rss)
                if rss>cfg["worker_rss_limit_mib"]*2^20; stopreason="external_rss_limit"; kill(proc); break; end
            end
            if time()-began>cfg["worker_case_timeout_seconds"]; stopreason="external_case_timeout"; kill(proc); break; end
        end
        wait(proc); close(logio)
        if isempty(stopreason)
            if !success(proc); print(read(logfile,String)); error("Worker failed; see $logfile"); end
            println("Completed DFS/BFS extension n=$n"); break
        end
        failure=Dict("run_id"=>"external_$(n)_$(current)","case_index"=>current,"nqubits"=>n,
          "algorithm"=>"worker","experiment"=>"external_failure","status"=>"resource_limited",
          "reason"=>stopreason,"peak_sampled_rss_mib"=>peak/2^20,"timestamp"=>string(now(UTC))*"Z")
        appendrow=open(outfile,"a"); println(appendrow,JSON3.write(failure)); close(appendrow)
        println("Stopped n=$n case $current: $stopreason; continuing"); skip=current
    end
end
include(joinpath(@__DIR__,"summarize_dfs_benchmarks.jl"))
