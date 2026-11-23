using JSON3, LinearAlgebra, Statistics, Random, Dates
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","FavorableCircuits.jl")); using .FavorableCircuits
include(joinpath(@__DIR__,"..","..","hierarchical_pauli","src","HierarchicalPauli.jl")); using .HierarchicalPauli
include(joinpath(@__DIR__,"..","..","hierarchical_pauli","src","HierarchicalPauliV2.jl")); using .HierarchicalPauliV2
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); outdir=joinpath(root,"results","favorable_v1"); mkpath(outdir); out=joinpath(outdir,"screen_v2.jsonl"); isfile(out)&&rm(out)
function emit(x); open(out,"a") do io; println(io,JSON3.write(x)); end end
function apply_map_sum(m,d)
    o=Dict{PauliWord,Float64}(); for (p,c) in d; for (q,v) in HierarchicalPauli.apply_map(m,p,c); o[q]=get(o,q,0.0)+v end end; HierarchicalPauli.merge_terms(o)
end
function cached_apply(c,p,b=4)
    d=Dict(p=>1.0); cache=Dict{UInt64,Any}(); built=0.0; chunks=0
    for last in length(c.gates):-b:1
        first=max(1,last-b+1); gs=c.gates[first:last]; key=UInt64(hash(gs)); m=get(cache,key,nothing)
        if m===nothing; z=@timed HierarchicalPauli.build_exact_block_map(gs); m=z.value; cache[key]=m; built+=z.time end
        d=apply_map_sum(m,d); chunks+=1
    end
    (;expectation=sum((v for (q,v) in d if q.x==0);init=0.0),terms=length(d),built_seconds=built,chunks,cache_size=length(cache))
end
function repeated_power(block,p,reps)
    m=HierarchicalPauli.build_exact_block_map(block); powers=[m]; pow=1
    while 2pow<=reps
        q=HierarchicalPauli.compose_maps(powers[end],powers[end]); sum((length(v) for v in values(q.columns));init=0)>5000 && break
        push!(powers,q); pow*=2
    end
    d=Dict(p=>1.0); used=0
    for i in length(powers):-1:1
        take=min(2^(i-1),reps-used); d=apply_map_sum(powers[i],d); used+=take; used==reps && break
    end
    (;expectation=sum((v for (q,v) in d if q.x==0);init=0.0),terms=length(d),powers=length(powers))
end
function run(f); GC.gc(); z=@timed f(); (;value=z.value,seconds=z.time,bytes=z.bytes) end
obs(c)=first(values(FavorableCircuits.observables(c.nqubits)))
families=(:echo,:repeated,:clifford_sandwich,:trotter,:periodic)
for family in families
  c=family==:echo ? family_case(family;n=4,depth=8,seed=11,perturb=:noncommuting,theta=.1) :
    family==:repeated ? family_case(family;n=3,depth=8,reps=8,seed=11) :
    family==:clifford_sandwich ? family_case(family;n=4,depth=8,reps=4,seed=11) :
    family==:trotter ? family_case(family;n=6,reps=4,theta=.05,perturb=.01) :
    family_case(family;n=3,reps=8,perturb=.01)
  p=obs(c); ref=upstream(c,p,:VectorPauliSum).expectation
  repeated_block=family==:repeated ? supported_local_block(3,8,11;kind=:mixed) : BestFirstPauli.SupportedGate[]
  methods=Dict{String,Function}(
    "B0_upstream"=>()->upstream(c,p,:VectorPauliSum),
    "B2_gatewise"=>()->HierarchicalPauli.propagate(c,p;block_size=1,threshold=1e-3),
    "B3_boundary"=>()->HierarchicalPauli.propagate(c,p;block_size=4,threshold=1e-3),
    "B4_cached_local"=>()->cached_apply(c,p,4),
    "B5_balanced_hier"=>()->HierarchicalPauliV2.v2_propagate(c.gates,p;strategy=:balanced,config=make_config(truncation_mode=:all_internal,coefficient_threshold=1e-3,max_support=4,max_columns=256)),
    "B6_support_hier"=>()->HierarchicalPauliV2.v2_propagate(c.gates,p;strategy=:support_constrained,config=make_config(truncation_mode=:all_internal,coefficient_threshold=1e-3,max_support=3,max_columns=64)),
    "B7_repeated_powers"=>()->(family==:repeated ? repeated_power(repeated_block,p,8) : (;status=:not_applicable,expectation=missing)),
    "B8_ondemand"=>()->HierarchicalPauliV2.v2_propagate(c.gates,p;strategy=:balanced,config=make_config(max_support=0,max_columns=1)))
  for (name,f) in methods
    run(f); for rep in 1:3
      z=run(f); r=z.value; ex=hasproperty(r,:expectation) ? r.expectation : missing
      emit((;family=String(family),method=name,rep,qubits=c.nqubits,gates=length(c.gates),reference=ref,expectation=ex,observed_error=ismissing(ex) ? missing : abs(ex-ref),seconds=z.seconds,allocated_bytes=z.bytes,terms=hasproperty(r,:terms) ? length(r.terms) : missing,status=hasproperty(r,:status) ? String(r.status) : "complete",timestamp=string(now(UTC))))
    end
  end
end
println(out)
