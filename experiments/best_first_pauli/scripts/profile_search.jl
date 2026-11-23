# Supplementary diagnostics: separate from benchmark timings and raw comparisons.
using Profile, JSON3
include(joinpath(@__DIR__,"..","src","BestFirstPauli.jl"))
using .BestFirstPauli
const ROOT = normpath(joinpath(@__DIR__,".."))
function workload(c,p)
    for _ in 1:20
        bestfirst(c,p;max_expansions=500_000,max_frontier_nodes=250_000,max_runtime_seconds=Inf)
    end
end
c=make_circuit(10,:ising,4,11); p=observable(10,:ZZ_center)
workload(c,p)
GC.gc(); Profile.clear(); Profile.init(n=10_000_000,delay=0.001)
Profile.@profile workload(c,p)
data,info=Profile.retrieve()
flat,lines=Profile.flatten(data,info)
open(joinpath(ROOT,"results","raw","profile_search.txt"),"w") do io
    Profile.print(IOContext(io,:displaysize=>(200,180)),flat,lines;format=:flat,sortedby=:count,mincount=5)
end
# A sample is one zero-delimited stack. Use metadata-free fetch to avoid treating
# thread/task/timestamp words as frames. Count each category at most once per stack.
stackdata=Profile.fetch(include_meta=false)
counts=Dict("samples"=>0,"heap_samples"=>0,"expand_samples"=>0,"search_samples"=>0)
stack=UInt[]
for ip in stackdata
    if ip==0
        isempty(stack) && continue
        frames=[frame for addr in stack for frame in get(info,addr,Base.StackTraces.StackFrame[])]
        counts["samples"]+=1
        any(f->occursin("frontier.jl",string(f.file)),frames) && (counts["heap_samples"]+=1)
        any(f->occursin("gates.jl",string(f.file)),frames) && (counts["expand_samples"]+=1)
        any(f->occursin("search.jl",string(f.file)),frames) && (counts["search_samples"]+=1)
        empty!(stack)
    else
        push!(stack,ip)
    end
end
write(joinpath(ROOT,"results","raw","profile_search.json"),JSON3.write((;
    circuit="n=10 ising layers=4 seed=11 ZZ_center",repetitions=20,max_expansions=500_000,
    sampling_interval_seconds=.001,counts,
    caveat="Stack-inclusive sample counts overlap and are noisy; separate instrumented run, not benchmark runtime.")))
println(counts)
