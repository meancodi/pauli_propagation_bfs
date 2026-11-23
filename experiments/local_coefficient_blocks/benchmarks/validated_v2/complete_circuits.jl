include(joinpath(@__DIR__,"common.jl"))
const OUT=joinpath(V2ROOT,"complete_circuits.jsonl")
open(OUT,"w") do io
    for n in (5,8), family in (:ising,:random), layers in (4,8), seed in (11,29)
        c=B.make_circuit(n,family,layers,seed); p=B.observable(n,:Z_center); data=Dict(p=>1.0)
        app=reverse(c.gates); up=upstream_objects(data,n,app)
        v=Float64[]; d=Float64[]; cross=0
        vt=PP.VectorPauliSum(up.pstrs); dt=from_terms(DenseLCB,data,[1,2])
        PP.propagate(up.ug,clone_upstream(vt),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
        x=deepcopy(dt); s=lcb_apply!(x,[[1,2]],[app]); cross=s.crossing
        for _ in 1:7
            push!(v,@elapsed z=PP.propagate(up.ug,clone_upstream(vt),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false))
            push!(d,@elapsed z=lcb_apply!(deepcopy(dt),[[1,2]],[app]))
        end
        rec=(;experiment="validated_complete_circuit",nqubits=n,family=String(family),layers,seed,
            primitive_gates=length(app),vector_seconds=median(v),dense_lcb_seconds=median(d),
            speedup=median(v)/median(d),crossing_gates=cross,final_expectation=expectation0(to_terms(x)),status="complete")
        JSON3.write(io,rec); write(io,'\n'); flush(io)
    end
end
