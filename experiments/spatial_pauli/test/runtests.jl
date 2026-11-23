using Test, Random, LinearAlgebra
include(joinpath(@__DIR__,"..","src","SpatialPauli.jl"))
using .SpatialPauli
const W=SpatialPauli.B

function compare(c, terms; sizes=(cld(c.nqubits,2),2))
    ref=sequential(c,terms)
    for size in sizes
        op=spatial(c,terms;leaf_size=size)
        got=expand_global(op)
        for p in union(keys(ref),keys(got))
            @test get(got,p,0.0) ≈ get(ref,p,0.0) atol=1e-10
        end
        @test expectation_zero(op) ≈ expectation_zero(ref) atol=1e-10
    end
end

@testset "CNOT orientation regression" begin
    for (a,b) in ((2,3),(3,2),(1,4),(4,1)), sa in (:I,:X,:Y,:Z), sb in (:I,:X,:Y,:Z)
        c=W.Circuit(4,[W.Clifford(:CNOT,a,b)])
        p=W.word([sa,sb],[a,b])
        ref=sequential(c,[p=>1.0])
        got=expand_global(spatial(c,[p=>1.0];leaf_size=2))
        for q in union(keys(ref),keys(got))
            @test get(got,q,0.0) ≈ get(ref,q,0.0) atol=1e-12
        end
    end
end

@testset "dense CNOT orientation" begin
    for (a,b) in ((2,3),(3,2),(1,4),(4,1)), sa in (:I,:X,:Y,:Z), sb in (:I,:X,:Y,:Z)
        g=W.Clifford(:CNOT,a,b); c=W.Circuit(4,[g]); p=W.word([sa,sb],[a,b])
        got=expand_global(spatial(c,[p=>1.0];leaf_size=2))
        u=W.dense_gate(g,4); exact=u' * W.dense_pauli(p,4) * u
        reconstructed=sum((v*W.dense_pauli(q,4) for (q,v) in got);init=zeros(ComplexF64,16,16))
        @test reconstructed ≈ exact atol=1e-12
    end
end

@testset "aligned spatial partitions" begin
    for n in (4,5,6,8,10,12,16), leaf_size in (1,2,3,4)
        root=partition(n,leaf_size); blocks=SpatialPauli.leaves(root)
        @test [(b.lo,b.hi) for b in blocks] == [(q,min(n,q+leaf_size-1)) for q in 1:leaf_size:n]
        @test vcat((collect(b.lo:b.hi) for b in blocks)...) == collect(1:n)
        @test root.lo==1 && root.hi==n
        for b in blocks
            @test SpatialPauli.lowest_common_ancestor(root,SpatialPauli.mask64(b)) == b
        end
        if length(blocks)>1
            s=(UInt64(1) << 0) | (UInt64(1) << (n-1))
            @test SpatialPauli.lowest_common_ancestor(root,s) == root
        end
    end
    @test [(b.lo,b.hi) for b in SpatialPauli.leaves(partition(12,2))] == [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12)]
    for n in (4,8,12)
        g=W.Rotation(W.word([:Z,:Z],[1,2]),0.4)
        op=spatial(W.Circuit(n,[g]),[W.word([:X],[1])=>1.0];leaf_size=2)
        @test op.cross_gates==0
    end
    @test SpatialPauli.mask64(partition(64,64))==typemax(UInt64)
    @test SpatialPauli.mask64(SpatialPauli.leaves(partition(64,2))[end]) == (UInt64(3) << 62)
end

@testset "cross rotations and phases" begin
    gens=[W.word([:X],[1]),W.word([:Y],[1]),W.word([:Z],[1]),
          W.word([:X,:X],[1,2]),W.word([:Y,:Y],[1,2]),W.word([:Z,:Z],[1,2]),
          W.word([:X,:Z],[1,2]),W.word([:X,:Y,:Z],[1,2,3])]
    syms=(:I,:X,:Y,:Z)
    for gen in gens, a in syms, b in syms, c in syms
        p=W.word([a,b,c],[1,2,3]); circ=W.Circuit(3,[W.Rotation(gen,0.37)])
        ref=sequential(circ,[p=>1.0]); got=expand_global(spatial(circ,[p=>1.0];leaf_size=1))
        for q in union(keys(ref),keys(got))
            @test get(got,q,0.0) ≈ get(ref,q,0.0) atol=1e-12
        end
    end
    for gen in gens[4:end], p in (W.word([:X],[1]),W.word([:Y],[2]),W.word([:Z],[3]))
        circ=W.Circuit(3,[W.Rotation(gen,0.37)])
        @test expectation_zero(spatial(circ,[p=>1.0];leaf_size=1)) ≈ W.dense_expectation(circ,p) atol=1e-12
    end
end

@testset "recursive spatial operator" begin
    for n in (4,5,6,8), seed in 1:3
        for family in (:ising,:random)
            c=W.make_circuit(n,family,1,seed)
            terms=[W.word([:Z],[cld(n,2)])=>1.0,W.word([:X,:X],[1,n])=>0.3]
            ref=sequential(c,terms)
            op=recursive_spatial(c,terms;leaf_size=2)
            got=recursive_expand(op)
            for p in union(keys(ref),keys(got))
                @test real(get(got,p,0im)) ≈ get(ref,p,0.0) atol=1e-10
                @test abs(imag(get(got,p,0im))) <= 1e-10
            end
            @test recursive_expectation(op) ≈ expectation_zero(ref) atol=1e-10
        end
    end
    for (a,b) in ((2,3),(3,2),(1,4),(4,1))
        c=W.Circuit(4,[W.Clifford(:CNOT,a,b)])
        p=W.word([:Y,:X],[a,b]); ref=sequential(c,[p=>1.0]); got=recursive_expand(recursive_spatial(c,[p=>1.0]))
        for q in union(keys(ref),keys(got)); @test get(got,q,0im) ≈ get(ref,q,0.0) atol=1e-12; end
    end
    c=W.Circuit(4,[W.Rotation(W.word([:X],[1]),0.37)])
    op=recursive_spatial(c,[W.word([:Z],[1])=>1.0])
    @test recursive_stats(op).reuse_events>0
    for n in (4,5,6,8)
        gates=W.SupportedGate[
            W.Rotation(W.word([:Z,:Z],[n,1]),0.31),
            W.Rotation(W.word([:X],[1]),0.27),
            W.Clifford(:CNOT,n,1),
            W.Rotation(W.word([:Y,:Y],[1,n]),0.19)]
        circ=W.Circuit(n,gates)
        terms=[W.word([:X,:Z],[1,n])=>0.7,W.word([:Z],[cld(n,2)])=>0.4]
        ref=sequential(circ,terms); got=recursive_expand(recursive_spatial(circ,terms;leaf_size=2))
        for p in union(keys(ref),keys(got))
            @test get(got,p,0im) ≈ get(ref,p,0.0) atol=1e-10
        end
    end
end

@testset "normalization and validation" begin
    a=SpatialPauli.Factor(W.word([:X],[1])=>1e-300,W.word([:Z],[1])=>-2.0)
    b=SpatialPauli.Factor(W.word([:I],[2])=>1.0)
    op=factorized(2,1,[Component(1.0,[a,b])])
    got=expand_global(op)
    @test got[W.word([:X],[1])] ≈ 1e-300 rtol=1e-12
    @test got[W.word([:Z],[1])] ≈ -2.0
    @test isempty(factorized(2,1,[Component(1.0,[a,b]),Component(-1.0,[a,b])]).components)
    @test_throws ArgumentError factorized(2,1,[Component(NaN,[a,b])])
    @test_throws ArgumentError factorized(2,1,[Component(1.0,[SpatialPauli.Factor(W.word([:X],[1])=>Inf),b])])
    @test_throws ArgumentError factorized(2,1,[Component(1.0,[SpatialPauli.Factor(W.word([:X],[1])=>1e-320,W.word([:Z],[1])=>1e100),b])])
    @test_throws ArgumentError SpatialPauli.propagate!(op,W.Clifford(:H,3))
    @test_throws ArgumentError SpatialPauli.propagate!(op,W.Rotation(W.word([:X],[3]),0.2))
    identity_circuit=W.Circuit(2,[W.Rotation(W.word([:I],[1]),0.31)])
    @test expand_global(spatial(identity_circuit,[W.word([:Z],[1])=>1.0];leaf_size=1)) == sequential(identity_circuit,[W.word([:Z],[1])=>1.0])
    @test recursive_expand(recursive_spatial(identity_circuit,[W.word([:Z],[1])=>1.0];leaf_size=1))[W.word([:Z],[1])] == 1.0
end

@testset "exact spatial propagation" begin
    I=W.word([:I]); X=W.word([:X],[1]); Z=W.word([:Z],[2])
    compare(W.Circuit(3,[]),[I=>2.0,X=>1.0,X=>-1.0,Z=>0.5];sizes=(1,2))
    a=SpatialPauli.Factor(W.word([:X],[1])=>1.0,W.word([:Z],[1])=>1.0)
    b=SpatialPauli.Factor(W.word([:X],[2])=>1.0,W.word([:Z],[2])=>1.0)
    op=factorized(2,1,[Component(1.0,[a,b])])
    @test length(op.components)==1
    @test length(expand_global(op))==4
    @test length(local_add(a,local_scale(-1.0,a)))==0
    compare(W.Circuit(4,[W.Rotation(W.word([:X],[1]),0.37)]),[W.word([:Z],[1])=>1.0];sizes=(2,))
    compare(W.Circuit(4,[W.Rotation(W.word([:Z,:Z],[1,2]),0.42)]),[W.word([:X],[1])=>1.0];sizes=(2,))
    compare(W.Circuit(4,[W.Rotation(W.word([:Z,:Z],[2,3]),0.42)]),[W.word([:X],[2])=>1.0];sizes=(2,))
    compare(W.Circuit(4,[W.Clifford(:CNOT,2,3)]),[W.word([:X,:Z],[2,3])=>1.0];sizes=(2,))
    gates=W.SupportedGate[W.Rotation(W.word([:X],[2]),0.21),W.Clifford(:CNOT,2,3),W.Rotation(W.word([:Z,:Z],[2,3]),0.31),W.Rotation(W.word([:Y],[3]),0.44)]
    compare(W.Circuit(4,gates),[W.word([:Z,:X],[2,3])=>1.0];sizes=(2,1))
    for n in (3,5,8), seed in 1:4
        c=W.make_circuit(n,:random,1,seed)
        compare(c,[W.word([:Z],[cld(n,2)])=>1.0];sizes=(cld(n,2),2))
    end
    for seed in 1:3
        c=W.make_circuit(4,:random,1,seed)
        p=W.word([:Z],[2])
        got=expectation_zero(spatial(c,[p=>1.0];leaf_size=2))
        @test got ≈ W.dense_expectation(c,p) atol=1e-10
        @test got ≈ W.upstream(c,p).expectation atol=1e-10
    end
    @test_throws ArgumentError SpatialPauli.propagate!(spatial(W.Circuit(2,[]),[I=>1.0]),:unsupported)
end
