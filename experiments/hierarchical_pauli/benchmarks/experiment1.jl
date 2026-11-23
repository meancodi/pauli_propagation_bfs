using JSON3, Statistics, Random, LinearAlgebra
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","HierarchicalPauli.jl")); using .HierarchicalPauli
LinearAlgebra.BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); out=joinpath(root,"results","experiment1.jsonl"); mkpath(dirname(out)); isfile(out)&&rm(out)
smoke=any(==("--smoke"),ARGS); lengths=smoke ? [2,4,8] : [2,4,8,16]
function emit(x); open(out,"a") do io; println(io,JSON3.write(x)); end end
for n in (2,3,4), len in lengths, seed in (11,29)
    c=make_circuit(n,:random,max(1,cld(len,2)),seed); gs=c.gates[1:min(len,length(c.gates))]
    for strategy in (:sequential,:left_deep,:balanced)
        GC.gc(); t=@timed (m=strategy==:sequential ? exact_sequential_map(gs) : recursive_map(gs;split=strategy); m)
        m=t.value; emit((;experiment=1,nqubits=n,block_length=length(gs),seed,strategy=String(strategy),support=length(m.support),generated=m.generated,unique_outputs=m.unique_outputs,cancelled=m.cancelled,entries=sum(length(v) for v in values(m.columns)),compression=m.generated/max(1,m.unique_outputs),seconds=t.time,allocated_bytes=t.bytes,package_version="PauliPropagation 0.9.2",julia_version=string(VERSION)))
    end
end
println(out)
