module LocalCoefficientBlocks

using Random
import Base: insert!
import PauliPropagation as PP
const B = let
    p = normpath(joinpath(@__DIR__, "..", "..", "best_first_pauli", "src", "BestFirstPauli.jl"))
    if !isdefined(Main, :BestFirstPauli)
        include(p)
        BestFirstPauli
    else
        Main.BestFirstPauli
    end
end

export DenseLCB, SparseLCB, PauliWord, outside_key, local_index,
       reconstruct_pauli, from_terms, from_upstream, to_terms, terms,
       insert!, remove_zeros!, propagate!, regroup!, nblocks, nterms,
       occupancy, support_mask, local_capacity

const PauliWord = B.PauliWord
const Term = Pair{PauliWord,Float64}

"Base-four local index: digits are I, X, Z, Y in active-qubit order."
@inline function local_digit(p::PauliWord, q::Int)
    b = UInt64(1) << (q - 1)
    x = !iszero(p.x & b); z = !iszero(p.z & b)
    (x ? 1 : 0) + (z ? 2 : 0)
end

@inline function support_mask(S::Vector{Int})
    m = UInt64(0)
    @inbounds for q in S
        m |= UInt64(1) << (q - 1)
    end
    m
end

@inline local_capacity(S) = 1 << (length(S) * 2)

@inline function outside_key(p::PauliWord, S::Vector{Int})
    m = support_mask(S)
    PauliWord(p.x & ~m, p.z & ~m)
end

@inline function local_index(p::PauliWord, S::Vector{Int})
    code = 0; stride = 1
    @inbounds for q in S
        code += local_digit(p, q) * stride
        stride <<= 2
    end
    code
end

function reconstruct_pauli(r::PauliWord, code::Integer, S::Vector{Int})
    x, z, t = r.x, r.z, Int(code)
    @inbounds for q in S
        d = t & 3; t >>= 2
        b = UInt64(1) << (q - 1)
        d & 1 != 0 && (x |= b)
        d & 2 != 0 && (z |= b)
    end
    PauliWord(x, z)
end

abstract type AbstractLCB end

mutable struct DenseLCB <: AbstractLCB
    active::Vector{Int}
    mask::UInt64
    blocks::Dict{PauliWord,Vector{Float64}}
end

mutable struct SparseLCB <: AbstractLCB
    active::Vector{Int}
    mask::UInt64
    blocks::Dict{PauliWord,Tuple{Vector{UInt16},Vector{Float64}}}
end

DenseLCB(S::AbstractVector{<:Integer}) = (a=sort(unique(Int.(S))); DenseLCB(a,support_mask(a),Dict{PauliWord,Vector{Float64}}()))
SparseLCB(S::AbstractVector{<:Integer}) = (a=sort(unique(Int.(S))); SparseLCB(a,support_mask(a),Dict{PauliWord,Tuple{Vector{UInt16},Vector{Float64}}}()))

function from_terms(::Type{DenseLCB}, data, S)
    out = DenseLCB(S)
    for (p,c) in data; insert!(out,p,c); end
    remove_zeros!(out)
end
function from_terms(::Type{SparseLCB}, data, S)
    out = SparseLCB(S)
    for (p,c) in data; insert!(out,p,c); end
    remove_zeros!(out)
end

"Convert any iterable of upstream `(paulistring, coefficient)` entries." 
function from_upstream(::Type{T}, upstream, S) where {T<:AbstractLCB}
    data = Pair{PauliWord,Float64}[]
    for item in upstream
        p, c = item
        # PauliPropagation's PauliString exposes symbols and sites.
        if p isa PauliWord
            push!(data, p => Float64(c))
        elseif p isa Integer
            n = hasproperty(upstream, :nqubits) ? getproperty(upstream, :nqubits) : maximum(S)
            str = PP.inttostring(p,n)
            push!(data, B.word(Symbol.(collect(str)),1:n) => Float64(c))
        else
            syms = hasproperty(p, :paulis) ? getproperty(p, :paulis) : getproperty(p, :string)
            sites = hasproperty(p, :sites) ? getproperty(p, :sites) : collect(1:length(syms))
            push!(data, B.word(Symbol.(syms), sites) => Float64(c))
        end
    end
    from_terms(T, data, S)
end

@inline nblocks(x::AbstractLCB) = length(x.blocks)
function nterms(x::DenseLCB)
    sum((count(!iszero, v) for v in values(x.blocks)); init=0)
end
function nterms(x::SparseLCB)
    sum((length(v[1]) for v in values(x.blocks)); init=0)
end
@inline occupancy(x::AbstractLCB) = nterms(x) / max(1, nblocks(x) * local_capacity(x.active))

function insert!(x::DenseLCB, p::PauliWord, c::Real)
    c == 0 && return x
    r = outside_key(p,x.active)
    v = get!(x.blocks,r) do
        zeros(Float64,local_capacity(x.active))
    end
    i = local_index(p,x.active) + 1; v[i] += Float64(c); x
end
function insert!(x::SparseLCB, p::PauliWord, c::Real)
    c == 0 && return x
    r = outside_key(p,x.active); pair = get(x.blocks,r,(UInt16[],Float64[])); inds, vals = pair
    i = UInt16(local_index(p,x.active)); j = searchsortedfirst(inds,i)
    if j <= length(inds) && inds[j] == i; vals[j] += Float64(c)
    else; insert!(inds,j,i); insert!(vals,j,Float64(c)); end
    x.blocks[r]=(inds,vals); x
end

function remove_zeros!(x::DenseLCB)
    for (r,v) in collect(x.blocks)
        all(iszero,v) && delete!(x.blocks,r)
    end
    x
end
function remove_zeros!(x::SparseLCB)
    for (r,(inds,vals)) in collect(x.blocks)
        keep = findall(!iszero,vals)
        x.blocks[r]=(inds[keep],vals[keep])
        isempty(keep) && delete!(x.blocks,r)
    end
    x
end

function terms(x::DenseLCB)
    out=Pair{PauliWord,Float64}[]
    for (r,v) in x.blocks, (i,c) in enumerate(v)
        c != 0 && push!(out,reconstruct_pauli(r,i-1,x.active)=>c)
    end
    out
end
function terms(x::SparseLCB)
    out=Pair{PauliWord,Float64}[]
    for (r,(inds,vals)) in x.blocks, (i,c) in zip(inds,vals)
        c != 0 && push!(out,reconstruct_pauli(r,i,x.active)=>c)
    end
    out
end
to_terms(x::AbstractLCB) = Dict(p=>c for (p,c) in terms(x))

function _local_gate(x::AbstractLCB,g)
    if g isa B.Rotation
        return (g.generator.x | g.generator.z) & ~x.mask == 0
    elseif g isa B.Clifford
        return g.kind != :CNOT ? (x.mask & (UInt64(1)<<(g.a-1)) != 0) :
            (x.mask & (UInt64(1)<<(g.a-1)) != 0 && x.mask & (UInt64(1)<<(g.b-1)) != 0)
    end
    false
end

function _prop_local!(x::DenseLCB,g)
    cap=local_capacity(x.active)
    for (r,v) in x.blocks
        out=zeros(Float64,cap)
        for (i,c) in enumerate(v)
            c==0 && continue
            for (q,cq) in B.expand(g,reconstruct_pauli(r,i-1,x.active),c)
                cq==0 && continue
                out[local_index(q,x.active)+1] += cq
            end
        end
        x.blocks[r]=out
    end
    remove_zeros!(x)
end
function _prop_local!(x::SparseLCB,g)
    for (r,(inds,vals)) in x.blocks
        outinds=UInt16[]; outvals=Float64[]
        for (i,c) in zip(inds,vals)
            for (q,cq) in B.expand(g,reconstruct_pauli(r,i,x.active),c)
                cq==0 && continue
                j=UInt16(local_index(q,x.active)); k=searchsortedfirst(outinds,j)
                if k<=length(outinds) && outinds[k]==j; outvals[k]+=cq
                else; insert!(outinds,k,j); insert!(outvals,k,cq); end
            end
        end
        x.blocks[r]=(outinds,outvals)
    end
    remove_zeros!(x)
end

"Apply an exact gate. Crossing gates use a correctness-first regrouping path." 
function propagate!(x::AbstractLCB,g)
    if _local_gate(x,g); _prop_local!(x,g)
    else
        data=to_terms(x); empty!(x.blocks)
        for (p,c) in data, (q,cq) in B.expand(g,p,c); cq != 0 && insert!(x,q,cq); end
        remove_zeros!(x)
    end
    x
end

function regroup!(x::Union{DenseLCB,SparseLCB}, S::AbstractVector{<:Integer})
    data=terms(x); a=sort(unique(Int.(S))); empty!(x.blocks)
    x.active=a; x.mask=support_mask(a)
    for (p,c) in data; insert!(x,p,c); end
    remove_zeros!(x)
end

end
