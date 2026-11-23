using Random, Statistics, JSON3
include(joinpath(@__DIR__,"..","..","src","LocalCoefficientBlocks.jl"))
using .LocalCoefficientBlocks
import PauliPropagation as PP
const B=LocalCoefficientBlocks.B
const V2ROOT=normpath(joinpath(@__DIR__,"..","..","results","validated_v2"))
mkpath(V2ROOT)
clone_upstream(x::PP.VectorPauliSum) = PP.VectorPauliSum(x.nqubits,copy(x.terms),copy(x.coeffs),x._terms_sorted)
clone_upstream(x::PP.PauliSum) = PP.PauliSum(x.nqubits,copy(x.terms))

safe_mask(n::Int) = n == 64 ? typemax(UInt64) : (UInt64(1) << n) - UInt64(1)

function validate_region(n,S)
    1 <= n <= 64 || throw(ArgumentError("n must be in 1:64"))
    all(q->1<=q<=n,S) || throw(ArgumentError("active qubit outside declared width"))
    length(unique(S))==length(S) || throw(ArgumentError("duplicate active qubit"))
end

function outside_from_code(code::UInt128,n::Int,S::Vector{Int})
    active=Set(S); x=UInt64(0); z=UInt64(0); t=code
    for q in 1:n
        q in active && continue
        d=Int(t & UInt128(3)); t >>= 2
        b=UInt64(1)<<(q-1)
        d&1!=0 && (x|=b); d&2!=0 && (z|=b)
    end
    PauliWord(x,z)
end

function synthetic_terms_v2(n::Int,S0,nt::Int,requested_occ::Float64,seed::Int)
    S=sort(Int.(S0)); validate_region(n,S)
    0 < requested_occ <= 1 || throw(ArgumentError("occupancy must be in (0,1]"))
    nt >= 0 || throw(ArgumentError("term count must be nonnegative"))
    cap=local_capacity(S); local_count=clamp(round(Int,cap*requested_occ),1,cap)
    maxgroups=BigInt(4)^(n-length(S))
    needed=cld(nt,local_count)
    BigInt(needed) <= maxgroups || throw(ArgumentError("infeasible: $needed outside groups exceeds $maxgroups"))
    rng=MersenneTwister(seed); codes=collect(UInt128(0):UInt128(needed-1)); shuffle!(rng,codes)
    out=Dict{PauliWord,Float64}(); per=Int[]
    for code in codes
        r=outside_from_code(code,n,S); order=collect(0:cap-1); shuffle!(rng,order)
        take=min(local_count,nt-length(out)); push!(per,take)
        for localidx in order[1:take]
            p=reconstruct_pauli(r,localidx,S)
            c=(rand(rng,Bool) ? 1.0 : -1.0)*(0.25+rand(rng))
            c*=rand(rng)<0.08 ? 1e-7 : 1.0
            out[p]=c
        end
    end
    mask=safe_mask(n)
    for p in keys(out)
        (p.x & ~mask)==0 && (p.z & ~mask)==0 || error("generator emitted out-of-range Pauli")
        outside_key(p,S).x & support_mask(S)==0 || error("active x bit leaked into outside key")
        outside_key(p,S).z & support_mask(S)==0 || error("active z bit leaked into outside key")
    end
    (;terms=out,n=n,active=S,capacity=cap,requested_terms=nt,actual_terms=length(out),
      requested_occupancy=requested_occ,actual_mean_occupancy=mean(per)/cap,
      occupancy_distribution=per ./ cap,outside_groups=length(per),seed)
end

function local_pattern(region::Vector{Int},g::Int)
    q1,q2=region
    base=B.SupportedGate[B.Clifford(:H,q1),B.Clifford(:CNOT,q1,q2),
        B.Rotation(B.word([:X,:Z],[q1,q2]),0.13)]
    [base[mod1(i,length(base))] for i in 1:g]
end

"Segments are in propagation order; the upstream circuit is reversed to match." 
function alternating_segments(g::Int,segments::Int)
    regs=[isodd(i) ? [1,2] : [3,4] for i in 1:segments]
    ss=[local_pattern(regs[i],g) for i in 1:segments]
    (regions=regs,segments=ss,chronological=reverse(vcat(ss...)))
end

function to_upstream(data,n)
    [PP.PauliString(n,B.symbols(p,n),collect(1:n),c) for (p,c) in data]
end

function from_upstream_terms(sum,n)
    d=Dict{PauliWord,Float64}()
    for (term,c) in sum
        s=PP.inttostring(term,n)
        p=B.word(Symbol.(collect(s)),1:n)
        d[p]=get(d,p,0.0)+Float64(c)
    end
    filter!(x->x.second!=0.0,d)
end

function upstream_objects(data,n,application_gates)
    pstrs=to_upstream(data,n)
    circuit=B.Circuit(n,reverse(application_gates))
    ug,theta,_=B.upstream_inputs(circuit,B.word([:Z],[1]))
    (;pstrs,ug,theta)
end

function lcb_apply!(x,regions,segments)
    local_time=0.0; regroup_time=0.0; crossing=0; changes=0
    for i in eachindex(segments)
        if x.active != regions[i]
            regroup_time += @elapsed regroup!(x,regions[i]); changes += 1
        end
        for g in segments[i]
            LocalCoefficientBlocks._local_gate(x,g) || (crossing += 1)
            local_time += @elapsed propagate!(x,g)
        end
    end
    (;local_time,regroup_time,crossing,changes)
end

function upstream_run(data,n,application_gates,T)
    pstrs=to_upstream(data,n)
    sum=T==:vector ? PP.VectorPauliSum(pstrs) : PP.PauliSum(pstrs)
    # The upstream adapter applies its gate list in reverse Heisenberg order.
    # Reverse the application-order list so both paths execute the same gates.
    circuit=B.Circuit(n,reverse(application_gates))
    ug,theta,_=B.upstream_inputs(circuit,B.word([:Z],[1]))
    PP.propagate(ug,sum,theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
end

function lcb_segments(T,data,regions,segments;count_crossing=true)
    x=from_terms(T,data,regions[1]); local_time=0.0; regroup_time=0.0; crossing=0; changes=0
    # `segments` are in application order, matching the initial region.
    for i in eachindex(segments)
        regions[i] == x.active || (regroup_time += @elapsed regroup!(x,regions[i]); changes += 1)
        for g in segments[i]
            LocalCoefficientBlocks._local_gate(x,g) || (crossing += 1)
            local_time += @elapsed propagate!(x,g)
        end
    end
    (;x,local_time,regroup_time,crossing,changes)
end

function expectation0(data)
    sum((c for (p,c) in data if p.x==0);init=0.0)
end
