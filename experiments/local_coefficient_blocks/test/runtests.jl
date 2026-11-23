using Test
include(joinpath(@__DIR__,"..","src","LocalCoefficientBlocks.jl"))
using .LocalCoefficientBlocks
const B = LocalCoefficientBlocks.B

function canonical(d)
    sort([(p.x,p.z,round(c,digits=12)) for (p,c) in d if abs(c)>1e-12])
end

@testset "LCB indexing and round trip" begin
    S=[4,1]
    for code in 0:15
        p=reconstruct_pauli(PauliWord(0,0),code,S)
        @test local_index(p,S)==code
        @test outside_key(p,S)==PauliWord(0,0)
    end
    data=Dict(B.word([:X,:Y,:Z],[1,4,6])=>0.5, B.word([:Z,:X],[2,5])=>-0.25)
    for T in (DenseLCB,SparseLCB)
        x=from_terms(T,data,S)
        @test canonical(to_terms(x))==canonical(data)
        @test nterms(x)==2
    end
end

@testset "exact local and crossing gates" begin
    data=Dict(B.word([:X,:Z],[1,4])=>0.7,B.word([:Y,:X],[2,5])=>-0.3)
    gs=B.SupportedGate[B.Clifford(:H,1),B.Clifford(:CNOT,1,4),B.Rotation(B.word([:X,:Z],[1,4]),0.37)]
    ref=copy(data)
    for g in gs
        out=Dict{PauliWord,Float64}()
        for (p,c) in ref, (q,cq) in B.expand(g,p,c)
            cq!=0 && (out[q]=get(out,q,0.0)+cq)
        end
        ref=Dict(p=>c for (p,c) in out if c!=0)
    end
    for T in (DenseLCB,SparseLCB)
        x=from_terms(T,data,[1,4])
        for g in gs; propagate!(x,g); end
        @test canonical(to_terms(x))==canonical(ref)
        regroup!(x,[2,5,6])
        @test canonical(to_terms(x))==canonical(ref)
    end
end

@testset "cancellation" begin
    p=B.word([:X],[1]); x=from_terms(DenseLCB,[p=>1.0,p=>-1.0],[1,2]); @test nterms(x)==0
    y=from_terms(SparseLCB,[p=>1.0,p=>-1.0],[1,2]); @test nterms(y)==0
end
