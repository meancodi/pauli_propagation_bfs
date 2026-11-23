using Random, Statistics, JSON3
include(joinpath(@__DIR__,"..","src","LocalCoefficientBlocks.jl"))
using .LocalCoefficientBlocks
const B=LocalCoefficientBlocks.B

const OUT=joinpath(@__DIR__,"..","results","synthetic.jsonl")
mkpath(dirname(OUT))

function synthetic_terms(n::Int, S::Vector{Int}, nt::Int, occ::Float64, seed::Int)
    rng=MersenneTwister(seed); cap=local_capacity(S)
    ng=max(1,ceil(Int,nt/max(1,round(Int,cap*occ))))
    groups=PauliWord[]; seen=Set{PauliWord}()
    while length(groups)<ng
        p=PauliWord(rand(rng,UInt64),rand(rng,UInt64))
        p=outside_key(p,S)
        p in seen || (push!(groups,p); push!(seen,p))
    end
    local_count=clamp(round(Int,cap*occ),1,cap)
    out=Dict{PauliWord,Float64}(); order=collect(0:cap-1)
    for r in groups
        shuffle!(rng,order)
        for code in order[1:local_count]
            length(out)>=nt && break
            p=reconstruct_pauli(r,code,S)
            # Mixed signs and a small near-canceling tail are intentional.
            c=(rand(rng)<0.5 ? -1.0 : 1.0)*(0.25+rand(rng))
            c*= rand(rng)<0.08 ? 1e-7 : 1.0
            out[p]=c
        end
        length(out)>=nt && break
    end
    out
end

function dict_step(data,g)
    out=Dict{PauliWord,Float64}()
    for (p,c) in data, (q,cq) in B.expand(g,p,c)
        cq!=0 && (out[q]=get(out,q,0.0)+cq)
    end
    filter!(x->x.second!=0.0,out)
end

function run_backend(data,gates,kind,S)
    x=kind==:dictionary ? copy(data) : from_terms(kind==:dense ? DenseLCB : SparseLCB,data,S)
    # Warm-up is deliberately outside the measured region.
    if kind==:dictionary; dict_step(x,gates[1])
    else; y=kind==:dense ? from_terms(DenseLCB,data,S) : from_terms(SparseLCB,data,S); propagate!(y,gates[1]); end
    GC.gc(); t=@elapsed begin
        for g in gates
            if kind==:dictionary; x=dict_step(x,g)
            else; propagate!(x,g); end
        end
    end
    alloc=@allocated begin
        for g in gates
            if kind==:dictionary; x=dict_step(x,g)
            else; propagate!(x,g); end
        end
    end
    (;seconds=t,allocated_bytes=alloc,final_terms=kind==:dictionary ? length(x) : nterms(x),
      blocks=kind==:dictionary ? missing : nblocks(x),
      final_occupancy=kind==:dictionary ? missing : occupancy(x))
end

function main()
    sizes=isempty(ARGS) ? [1_000,10_000,50_000,100_000] : parse.(Int,split(ARGS[1],","))
    occs=[0.05,0.15,0.30,0.60,0.90]; S=[1,2]; n=16; repeats=3
    gates_local=B.SupportedGate[B.Clifford(:H,1),B.Clifford(:S,2),B.Clifford(:CNOT,1,2),B.Rotation(B.word([:X,:Z],[1,2]),0.13)]
    gates_random=B.SupportedGate[B.Clifford(:H,3),B.Clifford(:CNOT,4,5),B.Rotation(B.word([:X,:Z],[6,7]),0.21),B.Clifford(:S,8)]
    open(OUT,"w") do io
        for nt in sizes, occ in occs
            if nt > 1_000_000; println(stderr,"skip $nt: resource cap"); continue; end
            data=synthetic_terms(n,S,nt,occ,1000+nt+round(Int,100*occ))
            actual=length(data); actual < nt && println(stderr,"requested $nt generated $actual")
            for pattern in (:local,:random), kind in (:dictionary,:dense,:sparse)
                gates=pattern==:local ? gates_local : gates_random
                vals=[run_backend(data,gates,kind,S) for _ in 1:repeats]
                rec=(;experiment="synthetic_fixed_region",requested_terms=nt,terms=actual,
                    requested_occupancy=occ,occupancy=actual/max(1,length(Set(outside_key.(keys(data),Ref(S))))*local_capacity(S)),
                    pattern=String(pattern),backend=String(kind),repetitions=repeats,
                    seconds=median(getfield.(vals,:seconds)),allocated_bytes=median(getfield.(vals,:allocated_bytes)),
                    final_terms=median(getfield.(vals,:final_terms)),julia_version=string(VERSION),threads=Threads.nthreads(),seed=1000+nt)
                JSON3.write(io,rec); write(io,'\n'); flush(io)
            end
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
