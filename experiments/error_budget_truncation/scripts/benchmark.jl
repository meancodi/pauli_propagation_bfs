using JSON3, Dates, LinearAlgebra, Statistics, Random
include(joinpath(@__DIR__,"..","..","best_first_pauli","src","BestFirstPauli.jl"))
using .BestFirstPauli
include(joinpath(@__DIR__,"..","src","ErrorBudgetTruncation.jl"))
using .ErrorBudgetTruncation
const B=BestFirstPauli
LinearAlgebra.BLAS.set_num_threads(1)

const ROOT=normpath(joinpath(@__DIR__,"..")); outdir=joinpath(ROOT,"results","raw")
mkpath(outdir); outfile=joinpath(outdir,"benchmark.jsonl")
append_results=any(==( "--append"),ARGS)
smoke=any(==("--smoke"),ARGS)
representative20=any(==( "--representative20"),ARGS)
ns=representative20 ? [20] : (smoke ? [5,10] : [5,10,20])
families=[:ising,:random,:cancellation]
layers=representative20 ? [1] : (smoke ? [1,2] : [1,2,4])
seeds=representative20 || smoke ? [11] : [11,29]
repeats=smoke || representative20 ? 2 : 3
budgets=[0.0,0.1,0.01,0.001]
policies=[:uniform,:front_loaded,:back_loaded,:cost_aware]
thresholds=[0.0,1e-3,1e-2,5e-2]

function cancellation_circuit(n,layers,seed)
    rng=MersenneTwister(seed); gs=B.SupportedGate[]
    for l in 1:layers
        for q in 1:n-1
            # Repeated near-opposite rotations deliberately create cancellation.
            push!(gs,Rotation(word([:Z,:Z],[q,q+1]),0.45))
            push!(gs,Rotation(word([:Z,:Z],[q,q+1]),-0.45+0.01*rand(rng)))
        end
        for q in 1:n
            push!(gs,Rotation(word([:X],[q]),0.32))
            push!(gs,Rotation(word([:X],[q]),-0.32+0.01*rand(rng)))
        end
    end
    Circuit(n,gs)
end

make(f,n,l,s)=f==:cancellation ? cancellation_circuit(n,l,s) : make_circuit(n,f,l,s)
function emit(row)
    open(outfile,"a") do io; println(io,JSON3.write(row)); end
end
!append_results && isfile(outfile) && rm(outfile)
for n in ns, family in families, layer in layers, seed in seeds
    c=make(family,n,layer,seed); obs=observable(n,:Z_center)
    exact = try upstream(c,obs,:PauliSum).expectation catch; missing end
    dense = n <= 8 ? dense_expectation(c,obs) : missing
    for budget in budgets, policy in policies
        propagate_budgeted(c,obs;epsilon_total=budget,policy=policy,max_seconds=60.0)
        for rep in 1:repeats
        GC.gc(); timed=@timed propagate_budgeted(c,obs;epsilon_total=budget,policy=policy,max_seconds=60.0)
        r=timed.value; observed=ismissing(exact) ? missing : abs(r.expectation-exact)
        emit((;algorithm="budgeted",nqubits=n,family=String(family),layers=layer,seed,observable="Z_center",
            policy=String(policy),epsilon_total=budget,min_abs_coeff=missing,repeat=rep,expectation=r.expectation,
            observed_error=observed,dense_error=ismissing(dense) ? missing : abs(r.expectation-dense),
            certified_bound=r.certified_bound,error_used=r.error_used,remaining_allowance=r.remaining_allowance,
            elapsed_ms=timed.time*1000,selection_ms=r.selection_seconds*1000,allocated_bytes=timed.bytes,
            gc_ms=timed.gctime*1000,peak_terms=r.peak_terms,retained_terms=r.retained_terms,
            discarded_terms=r.discarded_terms,status=String(r.status),reason=String(r.reason),
            package_version="PauliPropagation 0.9.2",julia_version=string(VERSION),threads=Threads.nthreads(),timestamp=string(now(UTC))))
        end
    end
    for threshold in thresholds
        fixed_threshold(c,obs;min_abs_coeff=threshold,max_seconds=60.0)
        for rep in 1:repeats
        GC.gc(); timed=@timed fixed_threshold(c,obs;min_abs_coeff=threshold,max_seconds=60.0)
        r=timed.value; observed=ismissing(exact) ? missing : abs(r.expectation-exact)
        emit((;algorithm="fixed_threshold",nqubits=n,family=String(family),layers=layer,seed,observable="Z_center",
            policy="fixed_threshold",epsilon_total=missing,min_abs_coeff=threshold,repeat=rep,expectation=r.expectation,
            observed_error=observed,dense_error=ismissing(dense) ? missing : abs(r.expectation-dense),
            certified_bound=missing,error_used=missing,remaining_allowance=missing,
            elapsed_ms=timed.time*1000,selection_ms=missing,allocated_bytes=timed.bytes,gc_ms=timed.gctime*1000,
            peak_terms=missing,retained_terms=r.retained_terms,discarded_terms=missing,status=String(r.status),reason="completed",
            package_version="PauliPropagation 0.9.2",julia_version=string(VERSION),threads=Threads.nthreads(),timestamp=string(now(UTC))))
        end
    end
    println("completed n=$n family=$family layer=$layer seed=$seed"); flush(stdout)
end
