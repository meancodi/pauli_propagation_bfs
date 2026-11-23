using Random, JSON3
include(joinpath(@__DIR__,"synthetic.jl")) # defines deterministic generator; main is not invoked on include
const OUT=joinpath(@__DIR__,"..","results","regrouping.jsonl"); mkpath(dirname(OUT))
open(OUT,"w") do io
    for nt in (1_000,10_000,50_000), gates_per_region in (1,2,4,8,16,32,64)
        data=synthetic_terms(16,[1,2],nt,0.3,7000+nt+gates_per_region)
        for T in (DenseLCB,SparseLCB)
            x=from_terms(T,data,[1,2]); S=[3,4]
            localg=B.SupportedGate[B.Clifford(:H,1),B.Clifford(:CNOT,1,2),B.Rotation(B.word([:X,:Z],[1,2]),0.11)]
            # Each segment has local work followed by a full regrouping.
            GC.gc(); t=@elapsed begin
                for _ in 1:8
                    for _ in 1:gates_per_region; for g in localg; propagate!(x,g); end; end
                    regroup!(x,S); S = S==[3,4] ? [1,2] : [3,4]
                end
            end
            rec=(;experiment="regrouping",terms=length(data),backend=string(nameof(T)),
                gates_per_region,segments=8,seconds=t,final_terms=nterms(x),julia_version=string(VERSION))
            JSON3.write(io,rec); write(io,'\n'); flush(io)
        end
    end
end
