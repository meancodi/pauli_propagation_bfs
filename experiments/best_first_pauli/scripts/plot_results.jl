using CSV, DataFrames, JSON3, Statistics, TOML
function consolidate_results()
    root = normpath(joinpath(@__DIR__,".."))
    raw = joinpath(root,"results","raw")
    rows = Dict{Symbol,Any}[]
    for f in sort(readdir(raw;join=true))
        endswith(f,".jsonl") && basename(f) != "failures.jsonl" || continue
        for line in eachline(f)
            isempty(strip(line)) && continue
            push!(rows,Dict(Symbol(k)=>(v===nothing ? missing : v) for (k,v) in JSON3.read(line)))
        end
    end
    isempty(rows) && error("No benchmark rows")
    df = DataFrame(rows)
    unique!(df,:run_id;keep=:last)
    df.reference = Vector{Union{Missing,Float64}}(missing,nrow(df))
    df.absolute_error = Vector{Union{Missing,Float64}}(missing,nrow(df))
    df.reference_source = fill("unavailable",nrow(df))
    keys=[:nqubits,:family,:layers,:seed,:observable]
    for group in groupby(df,keys)
        a=group[(group.backend .== "PauliSum") .& (group.status .== "exact"),:expectation]
        b=group[(group.backend .== "VectorPauliSum") .& (group.status .== "exact"),:expectation]
        dense=collect(skipmissing(group.dense_reference))
        ref=missing; source="unavailable"
        if !isempty(a) && !isempty(b)
            maximum(abs.(a .- first(b))) <= 1e-10 || error("Upstream backends disagree")
            ref=first(a); source="PauliSum + VectorPauliSum agreement"
        end
        if !isempty(dense)
            ismissing(ref) || abs(ref-first(dense)) <= 1e-10 || error("Dense reference disagrees")
            ref=first(dense); source="independent dense statevector + upstream"
        end
        for i in 1:nrow(group)
            group.reference[i]=ref
            group.absolute_error[i]=ismissing(ref) ? missing : abs(group.expectation[i]-ref)
            group.reference_source[i]=source
            if group.status[i] == "exact" && !ismissing(ref)
                group.absolute_error[i] <= 1e-10 || error("Exact result disagrees: $(group.run_id[i])")
            elseif group.backend[i] == "bestfirst" && !ismissing(ref)
                group.absolute_error[i] <= group.frontier_bound[i]+1e-10 || error("Frontier bound violated")
            end
        end
    end
    # External failures remain explicit in the primary CSV, even without an expectation.
    failures = joinpath(raw,"failures.jsonl")
    if isfile(failures)
        cfg = TOML.parsefile(joinpath(root,"benchmark_config.toml"))
        for line in eachline(failures)
            f=JSON3.read(line)
            cases=[(fam,l,s,o,0.0,"primary") for fam in cfg["families"] for l in cfg["layers"] for s in cfg["seeds"] for o in cfg["observables"]]
            f.backend == "bestfirst" && append!(cases,[(fam,4,11,"ZZ_center",b,"budget") for fam in cfg["families"] for b in cfg["budget_seconds"]])
            fam,l,s,o,b,suite=cases[f.case_index]
            row=Dict{Symbol,Any}(k=>missing for k in propertynames(df))
            merge!(row,Dict(:run_id=>"external_$(f.nqubits)_$(f.backend)_$(f.case_index)",:nqubits=>f.nqubits,
                :backend=>f.backend,:family=>fam,:layers=>l,:seed=>s,:observable=>o,:suite=>suite,
                :status=>f.status,:reason=>f.reason,:timestamp=>f.timestamp,:peak_rss_mib=>f.peak_sampled_rss_mib,
                :rss_source=>"external 50ms RSS samples until kill",:runtime_cap_seconds=>b))
            push!(df,row;promote=true)
        end
    end
    CSV.write(joinpath(raw,"benchmarks.csv"),df)
    stat(f,v)=isempty(collect(skipmissing(v))) ? missing : f(collect(skipmissing(v)))
    summary=combine(groupby(df,[:nqubits,:family,:layers,:seed,:observable,:backend,:suite,:runtime_cap_seconds]),
        nrow=>:repeats,:status=>(v->join(sort(unique(v)),"/"))=>:status,
        :elapsed_ms=>(v->stat(median,v))=>:median_ms,
        :elapsed_ms=>(v->stat(x->quantile(x,.75)-quantile(x,.25),v))=>:iqr_ms,
        :allocated_bytes=>(v->stat(median,v))=>:allocated_bytes,
        :peak_rss_mib=>(v->stat(maximum,v))=>:peak_rss_mib,
        :absolute_error=>(v->stat(median,v))=>:absolute_error,
        :frontier_bound=>(v->stat(median,v))=>:frontier_bound,
        :peak_frontier=>(v->stat(maximum,v))=>:peak_frontier,
        :peak_terms=>(v->stat(maximum,v))=>:peak_terms,
        :expanded_nodes=>(v->stat(median,v))=>:expanded_nodes,
        :completed_leaves=>(v->stat(median,v))=>:completed_leaves)
    CSV.write(joinpath(raw,"summary.csv"),summary)
    py=get(ENV,"PYTHON","python")
    run(`$py $(joinpath(@__DIR__,"plot_results.py"))`)
    println("Saved results/raw/benchmarks.csv, summary.csv, figures and report.md")
end
consolidate_results()
