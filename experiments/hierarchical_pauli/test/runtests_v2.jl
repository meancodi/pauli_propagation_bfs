using Test
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl")); using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","HierarchicalPauliV2.jl")); using .HierarchicalPauliV2

@testset "v2 exact strategies" begin
    c=make_circuit(3,:random,1,11); p=observable(3,:Z_center); ref=dense_expectation(c,p)
    for s in (:left_deep,:balanced,:support_constrained)
        r=v2_propagate(c.gates,p;strategy=s,config=make_config(truncation_mode=:none,max_support=3,max_columns=64))
        @test r.expectation ≈ ref atol=1e-10
        @test r.status == :complete
    end
end

@testset "genuine internal truncation and error bound" begin
    c=Circuit(1,[Rotation(word([:X]),0.7),Rotation(word([:Z]),0.8),Rotation(word([:Y]),0.9),Rotation(word([:X]),0.6)])
    p=word([:Z]); ref=dense_expectation(c,p)
    exact=v2_propagate(c.gates,p;config=make_config(truncation_mode=:none,max_support=1,max_columns=4))
    approx=v2_propagate(c.gates,p;config=make_config(truncation_mode=:all_internal,coefficient_threshold=.5,max_support=1,max_columns=4))
    @test approx.metrics.internal_truncations > 0
    @test approx.metrics.discarded_mass > 0
    @test approx.expectation != exact.expectation
    @test abs(approx.expectation-ref) <= approx.certified_bound + 1e-12
end

@testset "safe resource fallback and explicit cache" begin
    c=make_circuit(4,:random,1,29); p=observable(4,:Z_center)
    r=v2_propagate(c.gates,p;config=make_config(max_support=1,max_columns=1))
    @test r.status == :resource_limited
    @test r.metrics.resource_events > 0
    @test plan_summary(r.plan).kind == :sequential
    @test r.expectation ≈ dense_expectation(c,p) atol=1e-10
end
