# Best-first Pauli propagation

An independent CPU backend for fixed-angle Heisenberg expectation values on
`|0>^n`. The official registered PauliPropagation package is used unmodified for
the `PauliSum` and `VectorPauliSum` baselines. See [measured report](results/report.md).

## Setup and execution

Julia 1.10+ is required; the checked-in Manifest was resolved with Julia 1.13.4.
Use that version to reproduce this environment exactly. Python is used only for
standard Matplotlib figures. The Windows wrapper detects Juliaup's real executable
if the WindowsApps alias is inaccessible and automatically installs dependencies:

```powershell
.\experiments\best_first_pauli\scripts\setup.ps1
# Supply -Julia and -Python executable paths if automatic discovery fails.
# Use the Julia executable printed by setup in place of `julia` below.
julia --startup-file=no --threads=1 --project=experiments/best_first_pauli experiments/best_first_pauli/scripts/run_tests.jl
julia --startup-file=no --threads=1 --project=experiments/best_first_pauli experiments/best_first_pauli/scripts/run_benchmarks.jl
```

The setup wrapper sets `JULIA_DEPOT_PATH` to the repository's ignored `.julia-depot`
and `PYTHON` in the current shell. For a new shell, set those variables again or
rerun setup. Existing user depots can be reused by setting `JULIA_DEPOT_PATH`
before setup. No GPU is required. Do not run other CPU-intensive tasks during
measurements.

Cross-platform commands, from this experiment directory:

```sh
julia --startup-file=no scripts/setup.jl
python -m pip install -r requirements-plots.txt
julia --startup-file=no --threads=1 --project=. scripts/run_tests.jl
julia --startup-file=no --threads=1 --project=. scripts/run_benchmarks.jl 5
julia --startup-file=no --threads=1 --project=. scripts/run_benchmarks.jl 10
julia --startup-file=no --threads=1 --project=. scripts/run_benchmarks.jl 20
julia --startup-file=no --threads=1 --project=. scripts/plot_results.jl
julia --startup-file=no --threads=1 --project=. scripts/profile_search.jl
python scripts/validate_results.py
```

`PYTHON` can name a specific Python executable. Benchmark execution is offline;
run setup first. Setup checks the installed upstream source tree against the
Manifest's registered git-tree hash. Package provenance is saved in
`results/raw/upstream_verification.json`.

The benchmark script resumes completed configurations in existing per-backend
JSONL files. To perform a fresh run, move/archive `results/raw` first; do not mix
different configurations or source versions in one raw directory. Final CSV,
tables and plots are regenerated from the complete raw directory. Runtime/heap
caps and the grid are in [benchmark_config.toml](benchmark_config.toml).

## API

```julia
include("src/BestFirstPauli.jl")
using .BestFirstPauli
circuit = make_circuit(5, :ising, 2, 11)
result = bestfirst(circuit, observable(5, :Z_center);
    max_expansions=2_000_000, max_frontier_nodes=250_000,
    max_runtime_seconds=2.0)
@show result.expectation result.status result.reason result.frontier_bound

# Same physics and limits, deterministic LIFO traversal:
dfs_result = depthfirst(circuit, observable(5, :Z_center);
    max_expansions=2_000_000, max_frontier_nodes=250_000,
    max_runtime_seconds=2.0)
```

Supported gates are `Rotation(word([:X,:Z], [1,2]), theta)` and
`Clifford(:H,q)`, `Clifford(:S,q)`, `Clifford(:CNOT,control,target)`.
Build a `Circuit(n, gates)` in **Schrodinger order**. The backend traverses the
gates from last to first exactly once per path. The independent `(x,z)` encoding
uses qubit 1 in the least significant bit and supports up to 64 qubits; this does
not imply tractable computation at that size. Invalid/unsupported gates fail
explicitly. Observables are single Hermitian Pauli words with an optional signed
real `coefficient`; arbitrary complex non-Hermitian input operators are outside
this MVP. Complex multiplication phases are handled when deriving real branch
signs.

The max-heap orders paths by `abs(coeff)` and breaks ties by insertion sequence.
`priority_mode=:abs_coeff` is the default; a callable `coefficient -> score` can
replace ordering without changing the bound or terminal contributions. Scores
must be finite and nonnegative. No reference expectation enters prioritization.

`depthfirst` uses an explicit LIFO stack with the same `Path`, `expand` kernel,
terminal expectation rule, compensated sums, limits, and bound. The expansion
kernel returns cosine/unchanged before sine; DFS pushes in reverse so it visits
the cosine/unchanged child first. DFS has no merging, memoization, or pruning.
The algorithms therefore differ in frontier selection and associated data-
structure cost, rather than in the simulated physics.

Each nonterminal pop expands one gate and counts as one node expansion. Terminal
pops count as completed leaves. Zero floating-point children are omitted. No
coefficient cutoff is used; trigonometric roundoff at multiples of pi is retained.
Both the terminal sum and frontier mass use compensated accumulation. The reported
bound adds a small allowance for frontier bookkeeping roundoff; it is a numerical
bound with Float64 tolerance, not an interval-arithmetic proof.

- `:exact`: the frontier is empty, up to floating-point arithmetic.
- `:certified`: optional `error_tolerance` was satisfied using the frontier bound.
- `:resource_limited`: expansion, time, or frontier cap reached; `reason` identifies
  the cap and `frontier_bound` bounds omitted subtree contributions.

Before a branching expansion that would overflow memory, the algorithm leaves the
parent on the heap. Consequently no child is lost from the error ledger. Runtime
is checked between operations and may overrun by one operation. The path payload
is 48 bytes; 250,000 paths are about 11.44 MiB before vector spare capacity and
runtime overhead. The external process limit is deliberately separate from this
estimate.

## Reproducible circuits and comparisons

Open boundaries; initial state `|0>^n`; center index `cld(n,2)`; observables
`Z_center` and `Z_center Z_(center+1)`. Each Ising layer applies every nearest-neighbor
ZZ rotation, then X and Z rotations at each site, with angles respectively drawn
from `[0.31,0.35)`, `[0.47,0.51)`, `[0.23,0.27)` using `MersenneTwister(seed)`.
`R_Q(theta)=exp(-i theta Q/2)`. Random layers apply an H or S and a random
single-site Pauli rotation per site, followed by alternating disjoint nearest-
neighbor CNOT and random two-site Pauli rotations. Rotation angles lie in
`[-0.6,0.6)`. Seeds are 11 and 29. Explicit gate lists and exact parameter vectors
are saved, so future RNG changes do not erase the experiment specification.

Both upstream backends receive the same circuit via public `PauliRotation`,
`CliffordGate`, `PauliString` constructors and `propagate(...; thread=false,
min_abs_coeff=0.0, max_weight=Inf)`. All three use Float64 and one Julia/BLAS thread.
Conversion/setup are excluded from timing; propagation and final zero-state
overlap are included. One full configuration warmup precedes five separate
`@timed` runs with GC before each. Timing includes allocation and GC during the
run. Upstream peak terms are measured in a separate untimed `@peakpaulis` pass.

BFS searches unmerged paths while upstream returns a merged operator; path counts
and term counts are not equivalent. The report excludes interrupted runs from
speedup ratios. Budget experiments plot BFS error versus requested and measured
runtime; they do not claim matched-cost upstream comparisons.

The supervisor runs one backend/size worker at a time, samples RSS every 50 ms and
stops a process above 3 GiB or a case above 90 s including warmup/diagnostics. No
partial expectation is invented after an external kill. The failed case is logged
and later cases continue. Supported memory probes: Windows process-memory API,
Linux `/proc`, and macOS `ps`. Missing probes remain missing. Fast sampling cannot
guarantee prevention of every transient allocation spike; conservative internal
caps protect BFS independently.

## Files and metrics

- `src/`: encoding, kernels, heap, search, circuits, reference adapters, telemetry.
- `test/runtests.jl`: dense phase/sign tests, heap tests, ordering, upstream
  agreement, reproducibility, bounds and limits.
- `results/raw/benchmarks.csv`: repeat-level schema, missing cells for unavailable
  metrics; `summary.csv`: median/IQR tables; JSONL: original measurements.
- `results/raw/circuits/`: explicit circuit specifications; `metadata.json`:
  machine, threads, configuration and dependency hashes.
- `results/figures/`: standalone SVG and PNG plots; `results/report.md`: analysis.

RSS is process-lifetime high-water (including prior cases/JIT), allocation bytes
are total allocation traffic, and approximate node memory is live path payload.
These quantities must not be conflated. Heap occupancy median and isolated heap
time are unmeasured; mean occupancy and operation counts are recorded. A separate
warmed sampling profile (`profile_search.jl`) estimates heap stack prevalence
without affecting comparison timings. Tests use
`atol=1e-10` for circuit expectations and `1e-12` for elementary matrices.

Upstream API reference: [official repository](https://github.com/SparqleSim/PauliPropagation.jl).

## DFS versus best-first extension

The equal-budget extension is configured separately in
`dfs_benchmark_config.toml` and writes only under
`results_dfs_vs_bestfirst`, preserving every earlier result.

```powershell
# From the repository root, using the Julia executable printed by setup.ps1:
julia --startup-file=no --threads=1 --project=experiments/best_first_pauli experiments/best_first_pauli/scripts/run_dfs_benchmarks.jl 5
julia --startup-file=no --threads=1 --project=experiments/best_first_pauli experiments/best_first_pauli/scripts/run_dfs_benchmarks.jl 10
julia --startup-file=no --threads=1 --project=experiments/best_first_pauli experiments/best_first_pauli/scripts/run_dfs_benchmarks.jl 20
python experiments/best_first_pauli/scripts/validate_dfs_results.py
```

Omit the size arguments to run all three sizes. Runs resume after complete cases;
archive `results_dfs_vs_bestfirst` before changing source or configuration. The
grid uses expansion caps of 1k/10k/100k/1M, runtime caps of 1/10/100/1000 ms,
and certified frontier bounds of 1e-1 through 1e-4. Three post-warmup repeats are
recorded. Algorithm order alternates within each point. A separate instrumented
run measures kernel and push/pop bookkeeping; those timers are excluded from the
primary timing comparisons.

The supervisor imposes a 3 GiB sampled process-RSS guard and a 120-second per-case
guard. Search-level caps are 250,000 live paths, one million expansions, and one
second for accuracy-controlled runs. A killed or failed reference remains missing;
it is never inferred from a partial search. The primary research report is
[`results_dfs_vs_bestfirst/report.md`](results_dfs_vs_bestfirst/report.md), while
the earlier upstream comparison remains at [`results/report.md`](results/report.md).
