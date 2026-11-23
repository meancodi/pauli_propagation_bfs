using Test
include(joinpath(@__DIR__,"..","benchmarks","validated_v2","common.jl"))

@testset "bounded masked synthetic generator" begin
    for n in (4,8,16,32,64), S in ([1,2],[1,3])
        n < 3 && S==[1,3] && continue
        r=synthetic_terms_v2(n,S, min(8,BigInt(4)^(n-length(S))) |> Int, 0.15, n)
        mask=safe_mask(n)
        @test all((p.x & ~mask)==0 && (p.z & ~mask)==0 for p in keys(r.terms))
        @test all((outside_key(p,S).x & support_mask(S))==0 &&
                  (outside_key(p,S).z & support_mask(S))==0 for p in keys(r.terms))
        @test r.actual_terms==length(r.terms)
    end
    @test_throws ArgumentError synthetic_terms_v2(4,[1,2],17,0.05,1)
    @test_throws ArgumentError synthetic_terms_v2(4,[1,2],200,0.5,1)
end

@testset "matching alternating supports" begin
    d=synthetic_terms_v2(8,[1,2],100,0.6,7)
    a=alternating_segments(8,4)
    for (S,gs) in zip(a.regions,a.segments), g in gs
        @test LocalCoefficientBlocks._local_gate(from_terms(DenseLCB,d.terms,S),g)
    end
    x=lcb_segments(DenseLCB,d.terms,a.regions,a.segments)
    @test x.crossing==0
    @test x.changes==3
end

@testset "matched VectorPauliSum correctness" begin
    r=synthetic_terms_v2(8,[1,2],200,0.6,91); a=alternating_segments(4,2); app=vcat(a.segments...)
    up=upstream_objects(r.terms,8,app)
    v=PP.propagate(up.ug,PP.VectorPauliSum(up.pstrs),up.theta;min_abs_coeff=0.0,max_weight=Inf,thread=false)
    x=from_terms(DenseLCB,r.terms,a.regions[1]); lcb_apply!(x,a.regions,a.segments)
    @test expectation0(to_terms(x)) ≈ PP.overlapwithzero(v) atol=1e-10
    y=from_terms(SparseLCB,r.terms,a.regions[1]); lcb_apply!(y,a.regions,a.segments)
    @test expectation0(to_terms(y)) ≈ PP.overlapwithzero(v) atol=1e-10
end
