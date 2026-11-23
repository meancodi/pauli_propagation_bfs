include(joinpath(@__DIR__,"common.jl"))
const OUT=joinpath(V2ROOT,"realistic_intermediate.jsonl")

function capture(n,family,layers,seed)
    c=B.make_circuit(n,family,layers,seed); p=B.observable(n,:Z_center)
    cut=max(1,div(length(c.gates),2)); prefix=B.Circuit(n,c.gates[1:cut]); ug,theta,obs=B.upstream_inputs(prefix,p)
    initial=PP.PauliSum(obs)
    out=PP.propagate(ug,initial,theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
    from_upstream_terms(out,n),c,out
end

open(OUT,"w") do io
    for family in (:ising,:random), layers in (4,8), seed in (11,29)
        data,c,out=capture(8,family,layers,seed)
        app=alternating_segments(4,2); regions=app.regions; segs=app.segments; gates=vcat(segs...)
        up=upstream_objects(data,8,gates); vt=PP.VectorPauliSum(up.pstrs); dt=from_terms(DenseLCB,data,regions[1])
        PP.propagate(up.ug,clone_upstream(vt),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
        y=deepcopy(dt); lcb_apply!(y,regions,segs)
        vtimes=Float64[]; ltimes=Float64[]; conversion=Float64[]; ex=[]
        for _ in 1:5
            push!(vtimes,@elapsed z=PP.propagate(up.ug,clone_upstream(vt),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false))
            push!(ltimes,@elapsed z=lcb_apply!(deepcopy(dt),regions,segs))
            push!(conversion,@elapsed begin
                d=from_upstream_terms(PP.propagate(up.ug,clone_upstream(vt),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false),8)
                from_terms(DenseLCB,d,regions[1])
            end)
        end
        rec=(;experiment="captured_intermediate",family=String(family),layers,seed,
            captured_terms=length(data),captured_mean_occupancy=occupancy(from_terms(DenseLCB,data,[1,2])),
            continuation_gates=length(gates),vector_seconds=median(vtimes),dense_lcb_seconds=median(ltimes),
            capture_to_lcb_conversion_seconds=median(conversion),vector_expectation=PP.overlapwithzero(vt),
            status="complete")
        JSON3.write(io,rec); write(io,'\n'); flush(io)
    end
end
