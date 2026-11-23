using JSON3, Statistics, Printf
const ROOT=normpath(joinpath(@__DIR__,"..")); raw=joinpath(ROOT,"results","raw","benchmark.jsonl")
rows=[JSON3.read(line) for line in eachline(raw)]
function med(vals)
    good=[Float64(x) for x in vals if x !== nothing && x !== missing]
    isempty(good) ? missing : median(good)
end
open(joinpath(ROOT,"results","summary.md"),"w") do io
    println(io,"# Error-budgeted coefficient truncation\n")
    println(io,"Raw records: `results/raw/benchmark.jsonl`. Medians are over recorded repeats/cases; missing exact references are excluded.\n")
    println(io,"| algorithm | policy | n | family | layers | budget | median ms | median observed error | median certified bound | median selection % | incomplete |")
    println(io,"|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    keys=unique((String(r.algorithm),String(r.policy),Int(r.nqubits),String(r.family),Int(r.layers),ismissing(r.epsilon_total) ? "threshold" : string(r.epsilon_total)) for r in rows)
    for k in keys
        rs=filter(r->(String(r.algorithm),String(r.policy),Int(r.nqubits),String(r.family),Int(r.layers),ismissing(r.epsilon_total) ? "threshold" : string(r.epsilon_total))==k,rows)
        ms=med([r.elapsed_ms for r in rs]); err=med([r.observed_error for r in rs]); cert=med([r.certified_bound for r in rs])
        sel=med([r.selection_ms === nothing || r.selection_ms === missing ? missing : 100r.selection_ms/r.elapsed_ms for r in rs])
        inc=count(r->String(r.status)!="complete",rs)
        @printf(io,"| %s | %s | %d | %s | %d | %s | %.3f | %s | %s | %s | %d |\n",k...,ms,repr(err),repr(cert),repr(sel),inc)
    end
end
