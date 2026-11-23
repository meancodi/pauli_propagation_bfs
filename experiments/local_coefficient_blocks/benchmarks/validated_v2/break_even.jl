include(joinpath(@__DIR__,"common.jl"))

const OUT=joinpath(V2ROOT,"break_even.jsonl")
const REPEAT_SMALL=7
const REPEAT_LARGE=5

function median_mad(xs)
    m=median(xs); (;median=m,mad=median(abs.(xs .- m)),min=minimum(xs),max=maximum(xs))
end

function measure_case(meta, data, regions, segments, application_gates, repeats)
    n=meta.n
    up=upstream_objects(data,n,application_gates)
    templates=(vector=PP.VectorPauliSum(up.pstrs),pauli=PP.PauliSum(up.pstrs),
        dense=from_terms(DenseLCB,data,regions[1]),sparse=from_terms(SparseLCB,data,regions[1]))
    # Warm each exact path with an independent input copy.
    PP.propagate(up.ug,clone_upstream(templates.vector),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
    PP.propagate(up.ug,clone_upstream(templates.pauli),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
    y=deepcopy(templates.dense); lcb_apply!(y,regions,segments)
    y=deepcopy(templates.sparse); lcb_apply!(y,regions,segments)
    rows=NamedTuple[]
    for backend in (:VectorPauliSum,:PauliSum,:DenseLCB,:SparseLCB)
        steady=Float64[]; alloc=Int[]; localtimes=Float64[]; regroup=Float64[]; cross=Int[]; switches=Int[]
        e2e=Float64[]; outputs=Float64[]
        for _ in 1:repeats
            if backend==:VectorPauliSum || backend==:PauliSum
                t=@elapsed begin
                    out=PP.propagate(up.ug,clone_upstream(backend==:VectorPauliSum ? templates.vector : templates.pauli),up.theta;
                        min_abs_coeff=0.0,max_weight=Inf,thread=false)
                    e=PP.overlapwithzero(out)
                end
                push!(steady,t); push!(e2e,t); push!(outputs,e); push!(alloc,0); push!(localtimes,0.0); push!(regroup,0.0); push!(cross,0); push!(switches,0)
            else
                T=backend==:DenseLCB ? DenseLCB : SparseLCB
                x=deepcopy(backend==:DenseLCB ? templates.dense : templates.sparse)
                GC.gc()
                t=@elapsed stats=lcb_apply!(x,regions,segments)
                a=@allocated begin
                    z=deepcopy(backend==:DenseLCB ? templates.dense : templates.sparse)
                    lcb_apply!(z,regions,segments)
                end
                push!(steady,t); push!(e2e,t); push!(outputs,expectation0(to_terms(x))); push!(alloc,a)
                push!(localtimes,stats.local_time); push!(regroup,stats.regroup_time); push!(cross,stats.crossing); push!(switches,stats.changes)
            end
        end
        # Construction-inclusive timing is deliberately a separate contract.
        construction=Float64[]
        for _ in 1:repeats
            if backend==:VectorPauliSum || backend==:PauliSum
                t=@elapsed begin
                    ps=to_upstream(data,n); s=backend==:VectorPauliSum ? PP.VectorPauliSum(ps) : PP.PauliSum(ps)
                    out=PP.propagate(up.ug,s,up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false); PP.overlapwithzero(out)
                end
            else
                T=backend==:DenseLCB ? DenseLCB : SparseLCB
                t=@elapsed begin
                    x=from_terms(T,data,regions[1]); lcb_apply!(x,regions,segments); expectation0(to_terms(x))
                end
            end
            push!(construction,t)
        end
        push!(rows,(;meta...,backend=String(backend),repetitions=repeats,
            steady=median_mad(steady),end_to_end=median_mad(e2e),construction_inclusive=median_mad(construction),
            allocated_bytes=median(alloc),local_gate_seconds=median(localtimes),regroup_seconds=median(regroup),
            crossing_gates=maximum(cross),active_region_changes=maximum(switches),
            final_expectation=median(outputs),status="complete"))
    end
    rows
end

function main()
    sizes=isempty(ARGS) ? [1000,10000] : parse.(Int,split(ARGS[1],","))
    occs=parse.(Float64,split(get(ENV,"LCB_V2_OCC","0.15,0.60,0.90"),","))
    gs=parse.(Int,split(get(ENV,"LCB_V2_G","1,4,16,64"),","))
    segs=parse.(Int,split(get(ENV,"LCB_V2_SEGMENTS","2,8"),","))
    modes=Symbol.(split(get(ENV,"LCB_V2_MODES","zero_regroup,alternating"),","))
    n=16; seedbase=82000
    open(OUT,isfile(OUT) ? "a" : "w") do io
        for nt in sizes, occ in occs, g in gs, segments in segs, mode in modes
            S=[1,2]; r=try synthetic_terms_v2(n,S,nt,occ,seedbase+nt+g+segments) catch e
                rec=(;experiment="validated_v2",requested_terms=nt,requested_occupancy=occ,gates_per_segment=g,
                    segments,mode=String(mode),status="skipped",reason=sprint(showerror,e)); JSON3.write(io,rec); write(io,'\n'); continue
            end
            if mode==:zero_regroup
                regions=fill([1,2],segments); ss=[local_pattern([1,2],g) for _ in 1:segments]
            else
                a=alternating_segments(g,segments); regions=a.regions; ss=a.segments
            end
            app=vcat(ss...); repeats=nt<=10000 ? REPEAT_SMALL : REPEAT_LARGE
            meta=(experiment="validated_v2",requested_terms=nt,actual_terms=r.actual_terms,
                requested_occupancy=occ,actual_mean_occupancy=r.actual_mean_occupancy,
                outside_groups=r.outside_groups,active_qubits=length(S),n,
                gates_per_segment=g,primitive_gates=length(app),segments,mode=String(mode),seed=r.seed)
            for row in measure_case(meta,r.terms,regions,ss,app,repeats)
                JSON3.write(io,row); write(io,'\n'); flush(io)
            end
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
