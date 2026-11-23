using Test, Random, LinearAlgebra
include(joinpath(@__DIR__,"..","src","BestFirstPauli.jl"))
using .BestFirstPauli
const B = BestFirstPauli
LinearAlgebra.BLAS.set_num_threads(1)
@testset "Encoding and complex phases" begin
    words = [word([a,b]) for a in (:I,:X,:Y,:Z), b in (:I,:X,:Y,:Z)]
    for a in words, b in words
        @test word(symbols(a,2)) == a
        p,phase = product(a,b)
        am,bm = B.dense_pauli(a,2),B.dense_pauli(b,2)
        @test am*bm ≈ phase*B.dense_pauli(p,2) atol=1e-12
        @test commutes(a,b) == isapprox(am*bm,bm*am)
    end
    @test_throws ArgumentError word([:A])
    @test_throws ArgumentError word([:X,:Z],[1,1])
    @test word([:X],[64]).x == UInt64(1)<<63
end
@testset "Dense conjugations, rotations and Cliffords" begin
    gates = B.SupportedGate[Clifford(:H,1),Clifford(:S,1),Clifford(:CNOT,1,2),Clifford(:CNOT,2,1)]
    append!(gates,[Rotation(word([s],[q]),t) for s in (:X,:Y,:Z), q in 1:2, t in (0.0,pi/2,pi,-pi,0.31,-0.72)])
    for g in gates, a in (:I,:X,:Y,:Z), b in (:I,:X,:Y,:Z)
        p = word([a,b]); u = B.dense_gate(g,2)
        mat = sum(c*B.dense_pauli(q,2) for (q,c) in B.expand(g,p,1.0))
        @test mat ≈ u'*B.dense_pauli(p,2)*u atol=1e-12
        circ = Circuit(2,[g])
        for backend in (:PauliSum,:VectorPauliSum)
            @test bestfirst(circ,p).expectation ≈ upstream(circ,p,backend).expectation atol=1e-12
        end
    end
end
@testset "Heap order and stable ties" begin
    rng = MersenneTwister(4); h = B.Frontier()
    nodes = [B.Path(1,word([:Z]),1.0,Float64(rand(rng,1:8)),UInt64(i)) for i in 1:1000]
    foreach(p->B.heappush!(h,p),nodes)
    expected = sort(nodes;by=p->(-p.priority,p.seq_id))
    @test [B.heappop!(h) for _ in nodes] == expected
    @test isempty(h.nodes)
end
@testset "Order, exact agreement, seeds and resource bounds" begin
    c = Circuit(1,[Rotation(word([:Y]),0.43),Clifford(:H,1),Rotation(word([:X]),0.61)])
    p = word([:Z])
    @test bestfirst(c,p).expectation ≈ dense_expectation(c,p) atol=1e-12
    @test abs(dense_expectation(c,p)-dense_expectation(Circuit(1,reverse(c.gates)),p)) > 1e-3
    for family in (:ising,:random), seed in 1:6, layers in 1:2, obs in (:Z_center,:ZZ_center)
        c = make_circuit(5,family,layers,seed); p = observable(5,obs)
        @test c.gates == make_circuit(5,family,layers,seed).gates
        ref = dense_expectation(c,p)
        r = bestfirst(c,p)
        @test r.status == :exact
        @test r.frontier_bound == 0
        @test r.expectation ≈ ref atol=1e-10
        for backend in (:PauliSum,:VectorPauliSum)
            @test upstream(c,p,backend).expectation ≈ ref atol=1e-10
        end
        for budget in (0,1,10,100)
            r = bestfirst(c,p;max_expansions=budget)
            @test r.expanded_nodes <= budget
            @test abs(r.expectation-ref) <= r.frontier_bound+1e-10
        end
        r = bestfirst(c,p;max_frontier_nodes=2)
        @test r.peak_frontier <= 2
        @test abs(r.expectation-ref) <= r.frontier_bound+1e-10
        r2 = bestfirst(c,p;max_frontier_nodes=2)
        @test (r.expectation,r.expanded_nodes,r.frontier_bound) == (r2.expectation,r2.expanded_nodes,r2.frontier_bound)
    end
end
@testset "Boundary cases and validation" begin
    p=word([:Z]); c=Circuit(1,[Rotation(word([:X]),0.4)])
    @test bestfirst(Circuit(1,[]),p).expectation == 1
    @test bestfirst(c,p;coefficient=0).status == :exact
    @test bestfirst(c,p;max_runtime_seconds=0).reason == :runtime
    @test bestfirst(c,p;max_frontier_nodes=1).reason == :frontier
    @test bestfirst(c,p;max_expansions=0).frontier_bound >= 1
    @test bestfirst(c,p;error_tolerance=1.1).status == :certified
    @test bestfirst(c,p;error_tolerance=0).status == :exact
    certified=bestfirst(c,p;error_tolerance=0.4)
    @test certified.status == :certified
    @test certified.completed_leaves == 1
    @test abs(certified.expectation-dense_expectation(c,p)) <= certified.frontier_bound+1e-12
    @test bestfirst(c,p;priority_mode=x->1.0).expectation ≈ bestfirst(c,p).expectation
    @test bestfirst(Circuit(1,[Rotation(word([:X]),0.0)]),p).generated_children == 1
    @test bestfirst(Circuit(1,[Clifford(:H,1),Clifford(:S,1),Clifford(:H,1)]),p).completed_leaves == 1
    @test_throws ArgumentError Clifford(:T,1)
    @test_throws ArgumentError Circuit(1,["unknown"])
    @test_throws ArgumentError Circuit(1,[Clifford(:H,2)])
    @test_throws ArgumentError bestfirst(c,p;priority_mode=:unknown)
    @test_throws ArgumentError bestfirst(c,p;max_frontier_nodes=0)
    @test_throws ArgumentError bestfirst(c,p;priority_mode=x->NaN)
    @test_throws ArgumentError bestfirst(c,p;max_expansions=-1)
    @test_throws ArgumentError bestfirst(c,p;max_runtime_seconds=NaN)
    @test_throws ArgumentError bestfirst(c,p;error_tolerance=-1)
    @test_throws ArgumentError bestfirst(c,word([:Z],[2]))
    @test_throws ArgumentError Rotation(p,Inf)
end
@testset "DFS shared physics and deterministic traversal" begin
    for family in (:ising,:random), seed in 1:5, layers in 1:2, obs in (:Z_center,:ZZ_center)
        c=make_circuit(5,family,layers,seed); p=observable(5,obs)
        ref=upstream(c,p,:VectorPauliSum).expectation
        bfs=bestfirst(c,p); dfs=depthfirst(c,p)
        @test dfs.status == bfs.status == :exact
        @test dfs.expectation ≈ bfs.expectation atol=1e-10
        @test dfs.expectation ≈ ref atol=1e-10
        @test dfs.expanded_nodes == bfs.expanded_nodes
        @test dfs.generated_children == bfs.generated_children
        @test dfs.completed_leaves == bfs.completed_leaves
        @test dfs.frontier_bound == bfs.frontier_bound == 0
        for budget in (0,1,10,100)
            rb=bestfirst(c,p;max_expansions=budget)
            rd=depthfirst(c,p;max_expansions=budget)
            @test rb.expanded_nodes == rd.expanded_nodes
            @test rb.expanded_nodes <= budget
            if rb.status != :exact
                @test rb.expanded_nodes == budget
                @test rd.expanded_nodes == budget
            end
            @test abs(rb.expectation-ref) <= rb.frontier_bound+1e-10
            @test abs(rd.expectation-ref) <= rd.frontier_bound+1e-10
        end
        rd=depthfirst(c,p;max_frontier_nodes=2)
        @test rd.peak_frontier <= 2
        @test abs(rd.expectation-ref) <= rd.frontier_bound+1e-10
        rd2=depthfirst(c,p;max_frontier_nodes=2)
        @test (rd.expectation,rd.expanded_nodes,rd.frontier_bound,rd.completed_leaves) ==
              (rd2.expectation,rd2.expanded_nodes,rd2.frontier_bound,rd2.completed_leaves)
    end
    # Cosine/unchanged child is visited before sine: at the first leaf DFS has
    # followed the first child at every branch.
    c=Circuit(1,[Rotation(word([:X]),0.4),Rotation(word([:Y]),0.7)])
    r=depthfirst(c,word([:Z]);max_expansions=2)
    @test r.completed_leaves == 2
    @test r.stack_pushes == 5
    @test r.stack_pops == 4
    @test r.queue_pushes === missing
    @test r.frontier_nodes == 1
    @test r.expectation ≈ cos(0.7)*cos(0.4) atol=1e-12
    @test abs(r.expectation-dense_expectation(c,word([:Z]))) <= r.frontier_bound+1e-12
end
@testset "DFS limits, certification, and telemetry" begin
    c=make_circuit(5,:ising,2,11); p=observable(5,:ZZ_center); ref=dense_expectation(c,p)
    timed_dfs=depthfirst(c,p;max_runtime_seconds=0)
    timed_bfs=bestfirst(c,p;max_runtime_seconds=0)
    @test timed_dfs.reason == timed_bfs.reason == :runtime
    @test abs(timed_dfs.expectation-ref) <= timed_dfs.frontier_bound+1e-10
    @test abs(timed_bfs.expectation-ref) <= timed_bfs.frontier_bound+1e-10
    @test depthfirst(c,p;max_expansions=0).reason == :expansions
    @test depthfirst(c,p;max_frontier_nodes=1).reason == :frontier
    certified=depthfirst(c,p;error_tolerance=0.5)
    @test certified.status == :certified
    @test certified.frontier_bound <= 0.5
    @test abs(certified.expectation-ref) <= certified.frontier_bound+1e-10
    measured=depthfirst(c,p;max_expansions=100,measure_bookkeeping=true)
    @test measured.bookkeeping_seconds >= 0
    @test measured.kernel_seconds >= 0
    measuredb=bestfirst(c,p;max_expansions=100,measure_bookkeeping=true)
    @test measuredb.bookkeeping_seconds >= 0
    @test measuredb.kernel_seconds >= 0
    @test_throws ArgumentError depthfirst(c,p;max_expansions=-1)
    @test_throws ArgumentError depthfirst(c,p;max_frontier_nodes=0)
    @test_throws ArgumentError depthfirst(c,p;max_runtime_seconds=NaN)
    @test_throws ArgumentError depthfirst(c,p;error_tolerance=-1)
end
@testset "External RSS monitor" begin
    include(joinpath(@__DIR__,"..","src","telemetry.jl"))
    rss=process_rss(getpid())
    if Sys.iswindows() || Sys.islinux()
        @test !ismissing(rss)
        @test rss > 0
        @test rss <= Sys.maxrss()+2^20
    end
end
