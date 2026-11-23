using Test
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","FavorableCircuits.jl")); using .FavorableCircuits
@testset "five family generators" begin
    cases=[(:echo,echo_circuit(3,4,11;perturb=:noncommuting)),(:repeated,repeated_block_circuit(3,4,4,11)),(:clifford_sandwich,clifford_sandwich_circuit(3,4,2,11)),(:trotter,trotter_circuit(3,3,.05,.01)),(:periodic,periodic_circuit(3,4,.01))]
    for (name,c) in cases
        @test c.nqubits <= 3
        @test !isempty(c.gates)
        p=observable(c.nqubits,:Z_center)
        @test dense_expectation(c,p) ≈ upstream(c,p,:VectorPauliSum).expectation atol=1e-10
    end
end
@testset "inverse and periodic controls" begin
    c=echo_circuit(2,4,11;perturb=:identity); p=word([:Z],[1])
    @test isapprox(dense_expectation(c,p),1.0;atol=1e-10)
    pc=periodic_circuit(2,4,0.0;order=4); q=word([:X],[1])
    @test dense_expectation(pc,q) ≈ 0.0 atol=1e-10
end
