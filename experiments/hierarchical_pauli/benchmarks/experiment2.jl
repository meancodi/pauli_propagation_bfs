using JSON3, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","HierarchicalPauli.jl")); using .HierarchicalPauli
LinearAlgebra.BLAS.set_num_threads(1); root=normpath(joinpath(@__DIR__,"..")); out=joinpath(root,"results","experiment2.jsonl"); mkpath(dirname(out)); isfile(out)&&rm(out)
smoke=any(==("--smoke"),ARGS); ns=smoke ? [5,8] : [5,8,12,16]; bs=smoke ? [1,2,4] : [1,2,4,8,16]; ts=[1e-5,1e-3,1e-2]
function emit(x); open(out,"a") do io; println(io,JSON3.write(x)); end end
for n in ns, b in bs, threshold in ts
    c=make_circuit(n,:ising,2,11); p=observable(n,:Z_center); ref=n<=8 ? dense_expectation(c,p) : missing
    propagate(c,p;block_size=b,threshold=threshold); GC.gc(); z=@timed propagate(c,p;block_size=b,threshold=threshold); r=z.value
    emit((;experiment=2,nqubits=n,block_size=b,threshold,reference=ref,observed_error=ismissing(ref) ? missing : abs(r.expectation-ref),expectation=r.expectation,seconds=z.time,allocated_bytes=z.bytes,peak_terms=r.peak_terms,retained_terms=length(r.terms),truncation_seconds=r.truncation_seconds,truncation_events=r.truncation_events,dropped_mass=r.dropped_mass,status="complete"))
end
println(out)
