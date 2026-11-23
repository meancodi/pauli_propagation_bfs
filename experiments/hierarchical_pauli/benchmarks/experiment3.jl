using JSON3, LinearAlgebra
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","HierarchicalPauli.jl")); using .HierarchicalPauli
LinearAlgebra.BLAS.set_num_threads(1); root=normpath(joinpath(@__DIR__,"..")); out=joinpath(root,"results","experiment3.jsonl"); mkpath(dirname(out)); isfile(out)&&rm(out)
smoke=any(==("--smoke"),ARGS); ns=smoke ? [5,8] : [5,8,12]; strategies=(:sequential,:left_deep,:balanced,:support_constrained)
function emit(x); open(out,"a") do io; println(io,JSON3.write(x)); end end
for n in ns, strategy in strategies
    c=make_circuit(n,:random,2,29); p=observable(n,:Z_center); ref=n<=8 ? dense_expectation(c,p) : missing
    f=()-> strategy==:sequential ? propagate(c,p;block_size=4,threshold=1e-3) : propagate_map(c,p;strategy=strategy==:support_constrained ? :balanced : strategy,block_size=4,threshold=1e-3,max_support=smoke ? 4 : 6,max_columns=4096)
    f(); GC.gc(); z=@timed f(); r=z.value
    emit((;experiment=3,nqubits=n,strategy=String(strategy),reference=ref,observed_error=ismissing(ref) ? missing : abs(r.expectation-ref),expectation=r.expectation,seconds=z.time,allocated_bytes=z.bytes,peak_terms=r.peak_terms,map_count=get(r,:map_count,0),composition_seconds=get(r,:composition_seconds,0.0),dropped_mass=get(r,:dropped_mass,0.0),status="complete"))
end
println(out)
