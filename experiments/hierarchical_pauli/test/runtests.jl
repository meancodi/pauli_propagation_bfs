using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl"))
using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","HierarchicalPauli.jl"))
using .HierarchicalPauli

@testset "exact maps and ordering" begin
    c=make_circuit(4,:random,2,11); p=observable(4,:ZZ_center)
    seq=build_exact_block_map(c.gates)
    rec=recursive_map(c.gates;split=:balanced)
    left=recursive_map(c.gates;split=:left_deep)
    @test map_equal(seq,rec)
    @test map_equal(seq,left)
    out=apply_map_sum(seq,Dict(p=>1.0))
    up=upstream(c,p,:PauliSum)
    @test sum((v for (q,v) in out if q.x==0);init=0.0) ≈ up.expectation atol=1e-10
end

@testset "Clifford, phases, cancellation" begin
    c=Circuit(2,[Clifford(:H,1),Clifford(:S,1),Clifford(:CNOT,1,2)])
    m=build_exact_block_map(c.gates)
    @test map_equal(m,recursive_map(c.gates))
    x=word([:X],[1]); y=word([:Y],[1])
    @test merge_terms([(x,1.0),(x,-1.0),(y,0.25)]) == Dict(y=>0.25)
end

@testset "boundary truncation and references" begin
    c=make_circuit(5,:ising,2,29); p=observable(5,:Z_center)
    e=dense_expectation(c,p)
    for b in (1,2,4)
        r=propagate(c,p;block_size=b,threshold=0.0)
        @test r.expectation ≈ e atol=1e-10
    end
    @test propagate(c,p;block_size=2,threshold=0.2).truncation_events == cld(length(c.gates),2)
end
