using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl"))
using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","ErrorBudgetTruncation.jl"))
using .ErrorBudgetTruncation
const B=BestFirstPauli

@testset "deterministic certified selection" begin
    p0=B.PauliWord(0,0); px=B.word([:X]); py=B.word([:Y]); pz=B.word([:Z])
    s=Dict(px=>0.2,py=>-0.1,pz=>0.1,p0=>0.6)
    d=truncate_sum(s,0.1)
    @test d.discarded_mass == 0.1
    @test d.discarded_terms == 1
    @test d.retained_terms == 3
    @test haskey(d.retained,py) # tie resolved by (x,z): Y precedes Z, so Z is retained
    @test !haskey(d.retained,pz)
    @test truncate_sum(s,0.0).discarded_mass == 0.0
    @test truncate_sum(s,0.05).discarded_terms == 0
end

@testset "merge before truncation and ledger" begin
    c=B.Circuit(1,[B.Rotation(B.word([:X]),0.31),B.Rotation(B.word([:Y]),-0.22),B.Clifford(:H,1)])
    o=B.word([:Z])
    exact= B.upstream(c,o,:PauliSum).expectation
    r=propagate_budgeted(c,o;epsilon_total=0.0,policy=:uniform)
    @test r.status == :complete
    @test r.error_used == 0.0
    @test r.expectation ≈ exact atol=1e-12
    for policy in policy_names()
        x=propagate_budgeted(c,o;epsilon_total=0.01,policy=policy)
        @test x.error_used <= 0.01 + 16eps()
        @test x.error_used ≈ Base.sum(h.discarded_mass for h in x.history) atol=1e-15
        @test x.remaining_allowance ≈ 0.01-x.error_used atol=1e-14
    end
    # X + (-X) must disappear before a truncation decision.
    x=B.word([:X]); y=B.word([:Y]);
    @test ErrorBudgetTruncation.merged_sum([(x,1.0),(x,-1.0),(y,0.25)]) == Dict(y=>0.25)
end

@testset "budget policies, exact reference, and resource status" begin
    for policy in policy_names()
        @test isapprox(Base.sum(allocation_schedule(policy,7)),1.0;atol=1e-14)
        c=B.make_circuit(5,:ising,2,11); o=B.observable(5,:ZZ_center)
        r=propagate_budgeted(c,o;epsilon_total=0.0,policy=policy)
        @test r.expectation ≈ B.dense_expectation(c,o) atol=1e-10
    end
    c=B.make_circuit(5,:random,3,4); o=B.observable(5,:Z_center)
    limited=propagate_budgeted(c,o;epsilon_total=0.1,max_terms=1)
    @test limited.status == :incomplete
    @test limited.reason == :term_limit
end

@testset "fixed threshold is not a certified global budget" begin
    c=B.make_circuit(5,:ising,1,11); o=B.observable(5,:Z_center)
    r=fixed_threshold(c,o;min_abs_coeff=0.01)
    @test r.status == :complete
    @test r.total_discarded_mass >= 0
end
