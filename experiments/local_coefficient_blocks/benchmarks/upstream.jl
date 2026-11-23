using JSON3
include(joinpath(@__DIR__,"..","src","LocalCoefficientBlocks.jl"))
using .LocalCoefficientBlocks
const B=LocalCoefficientBlocks.B
const OUT=joinpath(@__DIR__,"..","results","upstream_circuits.jsonl"); mkpath(dirname(OUT))

function lcb_run(T,c,p,S)
    x=from_terms(T,Dict(p=>1.0),S)
    for g in reverse(c.gates); propagate!(x,g); end
    x
end
open(OUT,"w") do io
    for n in (5,8), family in (:ising,:random), depth in (4,8), seed in (11,29)
        c=B.make_circuit(n,family,depth,seed); p=B.observable(n,:Z_center); S=[1,2]
        # Warm-up each path once.
        B.upstream(c,p,:PauliSum); B.upstream(c,p,:VectorPauliSum)
        lcb_run(DenseLCB,c,p,S); lcb_run(SparseLCB,c,p,S)
        for (name, f) in (("PauliSum",()->B.upstream(c,p,:PauliSum)),
                          ("VectorPauliSum",()->B.upstream(c,p,:VectorPauliSum)),
                          ("DenseLCB",()->lcb_run(DenseLCB,c,p,S)),
                          ("SparseLCB",()->lcb_run(SparseLCB,c,p,S)))
            GC.gc(); t=@elapsed out=f(); rec=(;experiment="upstream_circuit",nqubits=n,
                family=String(family),depth=depth,seed=seed,backend=name,seconds=t,
                terms=name=="PauliSum" || name=="VectorPauliSum" ? out.unique_terms : nterms(out),
                expectation=name=="PauliSum" || name=="VectorPauliSum" ? out.expectation :
                    sum(c for (p,c) in terms(out) if p.x==0),package_version="PauliPropagation 0.9.2",julia_version=string(VERSION))
            JSON3.write(io,rec); write(io,'\n'); flush(io)
        end
    end
end
