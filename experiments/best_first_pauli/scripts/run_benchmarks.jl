using JSON3, TOML, Dates, SHA
include(joinpath(@__DIR__,"..","src","telemetry.jl"))
const ROOT = normpath(joinpath(@__DIR__,".."))
Threads.nthreads() == 1 || error("Run with --threads=1")
config = TOML.parsefile(joinpath(ROOT,"benchmark_config.toml"))
selected = isempty(ARGS) ? config["nqubits"] : parse.(Int,ARGS)
all(in(config["nqubits"]),selected) || error("Supported benchmark sizes: 5 10 20")
for p in ("raw","raw/circuits","work","figures")
    mkpath(joinpath(ROOT,"results",p))
end
metadata = Dict("julia_version"=>string(VERSION),"os"=>string(Sys.KERNEL),"machine"=>Sys.MACHINE,
    "cpu"=>Sys.cpu_info()[1].model,"logical_cpu_count"=>Sys.CPU_THREADS,"total_memory_bytes"=>Sys.total_memory(),
    "threads"=>Threads.nthreads(),"blas_threads"=>1,"timestamp"=>string(now(UTC))*"Z",
    "git_revision"=>try strip(read(`git -C $ROOT rev-parse HEAD`,String)) catch; "unknown" end,
    "gpu"=>"not queried; CPU-only; user reports NVIDIA RTX 4060 Laptop GPU",
    "config"=>config,"rss_poll_interval_seconds"=>0.05,
    "project_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"Project.toml")))),
    "manifest_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml")))))
metadata["source_sha256"] = Dict(relpath(joinpath(dir,file),ROOT)=>bytes2hex(sha256(read(joinpath(dir,file))))
    for (dir,_,files) in walkdir(joinpath(ROOT,"src")) for file in files if endswith(file,".jl"))
metadata["os_version"] = Sys.iswindows() ? try strip(read(`cmd /c ver`,String)) catch; "unavailable" end : try strip(read(`uname -a`,String)) catch; "unavailable" end
metadata["gpu"] = try strip(read(`nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader`,String)) catch; "not detected; CPU-only benchmark" end
oldmetadata = joinpath(ROOT,"results","raw","metadata.json")
if isfile(oldmetadata)
    old = JSON3.read(read(oldmetadata,String),Dict{String,Any})
    for field in ("config","manifest_sha256","source_sha256","julia_version")
        haskey(old,field) && old[field] != metadata[field] &&
            error("Existing results have different $field. Archive results/raw before starting a different experiment.")
    end
end
write(joinpath(ROOT,"results","raw","metadata_$(join(selected,'_')).json"),JSON3.write(metadata))
write(joinpath(ROOT,"results","raw","metadata.json"),JSON3.write(metadata))
for n in selected, backend in ("PauliSum","VectorPauliSum","bestfirst")
    outfile = joinpath(ROOT,"results","raw","$(n)_$(backend).jsonl")
    # Resume only fully recorded configurations; retain all earlier raw repeats.
    progressfile = joinpath(ROOT,"results","work","progress.json")
    skip = 0
    if isfile(outfile)
        previous = [JSON3.read(s) for s in eachline(outfile)]
        completed = [Int(r.case_index) for r in previous if r.repeat == config["repeats"]]
        skip = maximum(completed;init=0)
    end
    failurepath = joinpath(ROOT,"results","raw","failures.jsonl")
    if isfile(failurepath)
        failures = [JSON3.read(s) for s in eachline(failurepath)]
        skip = max(skip,maximum([Int(f.case_index) for f in failures if f.nqubits==n && f.backend==backend];init=0))
    end
    while true
        write(progressfile,JSON3.write((;case_index=skip+1,phase="startup"))*"\n")
        cmd = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$ROOT $(joinpath(@__DIR__,"benchmark_worker.jl")) $n $backend $skip $outfile $progressfile`
        logfile = joinpath(ROOT,"results","work","$(n)_$(backend)_after_$(skip).log")
        logio = open(logfile,"w")
        proc = run(pipeline(cmd;stdout=logio,stderr=logio);wait=false)
        worker_pid = getpid(proc)
        current = skip+1; began=time(); peakrss=0; reason=""
        while process_running(proc)
            sleep(0.05)
            state = try JSON3.read(last(readlines(progressfile))) catch; nothing end
            if state !== nothing && state.case_index != current
                current = state.case_index; began=time()
            end
            rss = process_rss(worker_pid)
            if !ismissing(rss)
                peakrss=max(peakrss,rss)
                if rss > config["worker_rss_limit_mib"]*2^20
                    reason="external_rss_limit"; kill(proc); break
                end
            end
            if time()-began > config["worker_case_timeout_seconds"]
                reason="external_case_timeout"; kill(proc); break
            end
        end
        wait(proc)
        close(logio)
        if isempty(reason)
            if !success(proc)
                print(read(logfile,String))
                error("Worker failed ($n $backend); see $logfile. No failed result is treated as exact.")
            end
            println("Completed n=$n backend=$backend; log: $logfile")
            break
        end
        failure=(;nqubits=n,backend,case_index=current,status=reason == "external_rss_limit" ? "resource_limited" : "timeout",
                  reason,peak_sampled_rss_mib=peakrss/2^20,timestamp=string(now(UTC))*"Z")
        open(joinpath(ROOT,"results","raw","failures.jsonl"),"a") do io
            println(io,JSON3.write(failure))
        end
        println("Stopped $n $backend case $current: $reason; continuing subsequent cases")
        skip=current
    end
end
include(joinpath(@__DIR__,"plot_results.jl"))
