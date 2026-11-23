# Design Document: Priority-Based Best-First Pauli Propagation

**Audience:** Codex coding agent  
**Status:** Implementation specification (MVP + benchmarks)  
**Language:** Julia  
**Primary comparison:** Official `PauliPropagation.jl` (`PauliSum` and `VectorPauliSum`)  
**Required system sizes:** 5, 10, and 20 qubits  
**Target development hardware:** Intel Core i7-13700HX, NVIDIA GeForce RTX 4060, 16 GB RAM , Windows laptop

## 1. Mission

Implement a **best-first path-exploration backend** for fixed-parameter Heisenberg-picture Pauli propagation. Instead of traversing paths depth-first or evolving a merged Pauli sum gate-by-gate, maintain a **max-priority frontier of partial paths** and expand the path with the highest priority first.

Benchmark it fairly against the **original, upstream PauliPropagation.jl** library on reproducible **5-, 10-, and 20-qubit** circuits, including exact-feasible and resource-bounded cases. Install and instantiate all project dependencies, run tests and benchmarks, and produce machine-readable results, plots, and a concise summary of when best-first wins or loses.

**Scope limitation:** This project initially evaluates fixed-angle expectation values, **not** complete LOWESA symbolic surrogate construction. Do not describe a fixed-angle heuristic as a parameter-uniform LOWESA guarantee.

## 2. Non-goals

- Do not replace or overwrite upstream PauliPropagation.jl internals.
- Do not introduce GPU acceleration, deep-learning-based priorities, symmetry-aware merging, or expensive preprocessing in the MVP. The target laptop has an NVIDIA RTX 4060, but CPU-only comparison is the required first stage; optional GPU experiments must be separate and justified by profiling.
- Do not assume best-first must outperform merging-based propagation; measure honestly.
- Do not silently use an incomplete truncated search as an exact result.
- Do not silently modify the upstream package or use a locally modified fork as the “original library” baseline.

## 3. Repository, environment, and installation

1. Inspect the existing repository and its project layout first. If working in an existing repository, add this as an isolated experiment (e.g. `experiments/best_first_pauli/`) without breaking existing functionality. If no repository is provided, create a self-contained Julia project.
2. Ensure Julia is installed and accessible (`julia --version`). Use a current compatible Julia version; upstream README specifies Julia **1.10+**. Do not assume package-specific internal APIs are stable.
3. Provide a **reproducible `Project.toml` and `Manifest.toml`**. Install via Julia `Pkg`; pin the resolved upstream revision/version in the manifest and record it in results. Prefer the registered upstream package `PauliPropagation` or the official upstream URL `https://github.com/SparqleSim/PauliPropagation.jl` on its default stable branch. Do not depend on a user fork for the official baseline.
4. Install only dependencies actually used. Likely dependencies: `PauliPropagation`, `DataStructures` (or an independently tested in-project binary max-heap), `BenchmarkTools`, `CSV`, `DataFrames`, `JSON3`, and a Julia plotting package if needed. Standard libraries `Random`, `Statistics`, `Dates`, `Test`, and `Pkg` need no external download.
5. Supply `scripts/setup.jl`, `scripts/test.jl` (or `test/runtests.jl`), and `scripts/run_benchmarks.jl` with documented command lines. Instantiate dependencies automatically and fail with actionable errors if installation fails. Do not request manual package installation for normal setup.
6. Record OS, CPU, Julia version, library version/revision, thread count, git revision, and run timestamp. Prioritize Windows (the expected target environment), while keeping Linux and macOS compatibility when practical. Verify the actual OS at runtime. Record the Intel Core i7-13700HX and NVIDIA RTX 4060 GPU if detected, but do not require GPU acceleration for these irregular search algorithms. Default to **single-threaded** runs for fair comparisons; optionally add separate multithreaded upstream results, clearly labeled.
7. No internet requests during benchmark execution beyond initial dependency installation.

**Upstream reference API (verify against installed version):**
```julia
using PauliPropagation
# circuit = ...; observable = ...; parameters = ...
psum = propagate(circuit, PauliSum(observable), parameters; thread=false)
value = overlapwithzero(psum)
# Also benchmark VectorPauliSum(observable) if supported by circuit/gates.
```
Consult current upstream examples for constructing Pauli rotations and benchmark circuits; do not guess missing constructors or private functions.

## 4. Mathematical scope and correctness contract

Let a circuit consist of unitary gates in **Schrödinger order** `G[1], ..., G[L]`; the Heisenberg propagation processes them **in reverse order**. Given an observable `O` and initial state `|0⟩^⊗n`, compute

\[
\langle O\rangle = \langle 0^{\otimes n}| U^\dagger O U |0^{\otimes n}\rangle.
\]

Each partial path contains an accumulated **signed real coefficient** `c`, its current Pauli word `P`, and its next gate index `k`. When the remaining gates have been processed, its contribution is:

\[
c\,\langle 0^{\otimes n}|P|0^{\otimes n}\rangle
=\begin{cases}c,&P\in\{I,Z\}^{\otimes n}\\0,&\text{otherwise.}\end{cases}
\]

For Pauli-rotation gates `R_Q(θ)=exp(-i θ Q/2)`:

- When `[P,Q]=0`, produce **one unchanged child**.
- When `{P,Q}=0`, produce the correct cosine and sine branches with the **correct real phase/sign** under the library's Heisenberg convention. Do **not** hardcode signs from an unverified convention; validate elementary cases against dense matrix calculations and upstream.
- Clifford gates produce a single signed Pauli child.
- Skip mathematically zero-coefficient children.
- **Initial supported gate set:** the intersection of upstream gates needed for the chosen circuits (Pauli rotations and supported Clifford gates). If a gate cannot be implemented accurately, fail clearly instead of silently dropping it.

For 5/10/20 qubits, use an efficient Pauli encoding (e.g. `(xmask::UInt64, zmask::UInt64)` or the package's documented stable representation). Include phase/sign handling and correctness tests; do not confuse the user's old dibit encoding with the symplectic `(x,z)` representation.

## 5. Best-first algorithm

### 5.1 Path representation

Store each path with:

- `next_gate`: index into reversed gate sequence;
- `pauli`: compact encoded Pauli string;
- `coeff`: signed coefficient (default `Float64`);
- `priority`: `Float64`;
- `seq_id`: monotonic insertion counter for deterministic tie-breaking.

Avoid storing complete path histories in memory. Optionally store lightweight counters for diagnostics. Keep per-node storage compact and type-stable.

### 5.2 Priority function

**MVP default:** `priority = abs(coeff)`.

Interpret this as an upper bound on the absolute **net final expectation contribution of that unfinished subtree** for fixed gate parameters, normalized initial state, and unitary remaining evolution:

\[
|c\,\mathrm{Tr}(\rho V^\dagger PV)| \le |c|.
\]

This is a safe bound but may be very loose. Do not claim it predicts actual contribution. Add an extensible `priority_mode` interface (initially only `:abs_coeff`; optional `:abs_coeff_weight` clearly marked heuristic). No machine-learning model in MVP.

**Critical invariant:** Path priorities determine expansion order; they do not change the exact result if all paths are explored.

### 5.3 Expansion loop (pseudocode)

```text
frontier = maxheap([Path(next_gate=L, pauli=O, coeff=1.0)])
completed_sum = 0
completed_leaves = 0
expanded_nodes = 0
while frontier not empty:
    if resource limit reached: break with status=:incomplete
    p = pop_max(frontier)
    if p.next_gate == 0:
        completed_sum += p.coeff * overlap_zero(p.pauli)
        completed_leaves += 1
        continue
    expanded_nodes += 1
    gate = reversed_gates[p.next_gate]
    children = heisenberg_expand(gate, p.pauli, p.coeff)
    for child in children:
        if child.coeff != 0:
            push frontier(child with next_gate=p.next_gate-1,
                          priority=abs(child.coeff))
    update telemetry and resource checks
end
return completed_sum, frontier_bound, stats, status
```

Implement indexing carefully and test with 1–3-gate hand-built circuits. Choose **node expansions** as the primary budget unit, not only completed leaves.

### 5.4 Correctness of interruption and optional accuracy control

At any stopping point, the unfinished frontier partitions the omitted search into disjoint subtrees. A safe **absolute error bound** is:

\[
E_{\mathrm{frontier}}=\sum_{p\in frontier}|c_p|.
\]

Therefore `abs(exact_expectation - completed_sum) <= E_frontier` in exact arithmetic, assuming no unexplained pruning, correct propagation, and exact arithmetic; in Float64, allow a numerical tolerance. Track `frontier_abs_coeff_sum` incrementally (add on push, subtract on pop) to avoid scanning the heap each step.

- `status=:exact` only when the frontier is empty.
- `status=:certified` only if an **optional** `error_tolerance` was met using the valid frontier bound.
- `status=:resource_limited` when budget or memory limits are hit first; report the remaining bound, do not claim requested accuracy was met.
- A low priority score is **not** permission to delete a path without including it in a persistent discarded-subtree error ledger.
- Do not introduce truncation in the main exact-search comparison; add a separately labeled approximate experiment if useful.

### 5.5 Memory controls for the target 16 GB RAM laptop

Support `max_frontier_nodes`, `max_expansions`, and optional `max_runtime_seconds`; choose conservative defaults and make the 20-qubit cases safely interruptible.

When `max_frontier_nodes` would be exceeded, **stop and report incomplete** in MVP (rather than silently discarding paths). Optionally implement a later memory-bounded hybrid DFS continuation, but do not make that a prerequisite.

Display the current frontier, peak frontier, approximate node memory, RSS where obtainable, and whether the run ended due to resource limits. Provide a conservative **16 GB RAM** default profile, leaving substantial RAM for the OS, Julia runtime, upstream comparison processes, and other applications. Make caps configurable, and avoid OS swapping. Don't pretend RSS is the same as Julia heap allocation.

## 6. Baselines: original PauliPropagation.jl

Compare best-first against **both** standard upstream backends where supported:

1. **`PauliSum`** via the public `propagate` API.
2. **`VectorPauliSum`** via the same public API; this may be the stronger official CPU baseline.

Match exactly:

- Circuit structure, gate sequence and angle values.
- Initial observable, qubit ordering, and initial `|0⟩^⊗n` state.
- Float precision and supported gate set.
- Single thread for primary results (`thread=false` if supported).
- No truncation for **exact-feasible** cases; use the same stopping/accuracy criteria only when a genuinely comparable mode exists.

**Fairness caveat:** Upstream `propagate()` is a merged, gate-by-gate algorithm. It does **not** use best-first's partial-tree frontier budget. Therefore, `max_expansions`/`max_frontier_nodes` are **not equivalent** to upstream's term count or `max_weight`. Report such resource-limited runs separately; do not label them as matched-cost head-to-head comparisons.

For bounded approximate comparisons, use **equal runtime and/or equal peak-memory budgets**, compare achieved error versus a common reference, and explicitly state the baseline's truncation settings (e.g. `min_abs_coeff`, `max_weight`). Do not force an upstream run to finish if memory pressure is unsafe.

As an optional supplementary reference, implement a tiny unmerged DFS enumerator sharing the same expansion kernel to separate traversal-order overhead from the impact of merging. This is useful because the official upstream baseline is primarily merging-based, not raw DFS.

## 7. Benchmark circuits: 5, 10, 20 qubits

Create reproducible circuits from **documented upstream constructors** where possible. Include at least two families:

**Family A — 1D tilted-field Ising Trotter circuits:** local `ZZ` nearest-neighbor interactions and noncommuting single-qubit `X`/`Z` rotations; choose deterministic angles and an observable such as `Z` on the center qubit. The user-facing size levels are exactly `n ∈ {5,10,20}`. Record boundary conditions and the precise gate/Trotter conventions.

**Family B — Random local Clifford + Pauli-rotation circuits:** seeded random local Clifford layers and parameterized single/two-qubit Pauli rotations, chosen to produce varied branching and reconvergence. Avoid unsupported gates.

For each `n`, run a small pilot grid of depth/layer settings, e.g. 1, 2, and 4 Trotter layers, **but** adapt down if raw path counts or resource estimates become intractable. Preserve at least one runnable case per family and qubit count. Never silently skip 20 qubits; record `resource_limited`, `timeout`, or `not_feasible`, with reason. Avoid presenting 20-qubit results as exact if the full computation did not finish.

Use fixed parameter vectors and seeds saved to results. Prefer `|0⟩^⊗n`, `Z_center`, and at least one additional observable to ensure the evaluation is not accidentally dominated by zero contributions.

## 8. Correctness testing

Mandatory tests:

1. Pauli encoding/decoding, identity, phase signs, commutation/anticommutation for `I,X,Y,Z` and multi-qubit strings.
2. Single-qubit rotations (`Rx`, `Ry`, `Rz`) and Clifford (`H`, `S`, `CNOT`) against analytic/dense-matrix conjugations.
3. Reverse-gate ordering: 2–3 noncommuting gates where incorrect order changes the result.
4. Exact best-first versus upstream exact `PauliSum` and `VectorPauliSum` on tractable random circuits.
5. For 5-qubit small circuits, independently cross-check expectation using a dense state-vector or direct matrix computation (no need for large external dependency).
6. Incomplete-run bound: verify `|reference - partial_result| <= frontier_bound + numerical_tolerance` on many small randomized cases.
7. Tie-handling reproducibility; empty frontier; zero children; zero/π/π/2 rotation angles; all-Clifford circuits; unsupported gate error.
8. Ensure no leakage of expected result from reference into best-first prioritization.

Use explicit numerical tolerances (e.g. `atol=1e-10` for small tractable Float64 circuits, adjusted only with justification). Never remove failing tests to get green CI.

## 9. Metrics and data collection

**Primary:**

- Wall-clock simulation time (exclude Julia compilation, setup, plotting, and reference calculation).
- Accuracy: `abs(value - reference)`, plus frontier error bound when relevant.
- Peak resident memory (RSS) if available; label source and sampling method.
- Peak frontier nodes (best-first) / peak term count (upstream); do not treat these as identical measures.

**Secondary:**

- Total node expansions, generated children, completed leaves, queue pushes/pops, maximal tree depth.
- Peak/median heap size; total heap-operation time or sampled overhead.
- Julia allocated bytes and GC time.
- Unique Pauli terms at completion, where available.
- Status, runtime cap, node cap, and any timeout.
- Upstream `@countpaulis` or `@peakpaulis` if supported by the resolved version.

**Data schema:** one row per repeat, with `run_id`, `algorithm`, `backend`, `nqubits`, `family`, `layers`, `seed`, `observable`, `theta_spec`, `truncation_config`, `priority_mode`, `status`, `expectation`, `reference`, `absolute_error`, `frontier_bound`, `elapsed_ms`, `allocated_bytes`, `gc_ms`, `peak_rss_mib`, `peak_frontier`, `peak_terms`, `expanded_nodes`, `generated_children`, `completed_leaves`, `julia_version`, `package_revision`, `threads`, and `timestamp`. Unavailable values should be **missing**, never invented zeros.

## 10. Benchmark procedure and presentation

1. **Warm up** each benchmark configuration before timing to exclude Julia JIT compilation from measured runs.
2. Use at least **5 measured repetitions** of feasible cases and report median and spread (IQR or min/max); resource-heavy 20-qubit cases may use fewer repetitions, explicitly noted.
3. Run a **single-threaded primary suite** on the same machine. Optional separate multithreaded suite is allowed but never mix it into single-thread comparisons.
4. Produce tables for each `n=5,10,20`, each circuit family/depth, and every algorithm/backend.
5. Generate charts: (a) median runtime by algorithm, (b) peak RSS/allocated memory by algorithm, (c) absolute error versus wall-clock budget for budgeted runs, (d) peak frontier/term growth, and (e) performance ratio `upstream_time / best_first_time` **only when comparing fully completed equivalent tasks**.
6. Mark **incomplete** or **resource-limited** cases clearly rather than excluding them from the report. Avoid misleading speedup claims if one algorithm is approximate while another is exact.
7. Provide a concise report covering cases in which best-first wins, loses, and runs out of memory, including priority-queue overhead and implications for 16 GB RAM.

## 11. Implementation milestones

### M1 — Environment and reference baseline
- Setup Julia environment and dependencies.
- Create deterministic circuits and run official upstream `PauliSum`/`VectorPauliSum` on 5-qubit smoke tests.
- Record upstream version/revision and confirm exact public API.

### M2 — Best-first correctness
- Implement compact path representation, max-heap, Heisenberg expansion kernels, and `overlapwithzero`-equivalent terminal overlap.
- Validate elementary gates, ordering, and exact expectation values.
- Implement tests for `frontier_abs_coeff_sum` and interruption.

### M3 — Resource tracking and limits
- Add node/time/frontier caps and telemetry.
- Ensure graceful resource-limited exits, especially for 20 qubits on 16 GB RAM.

### M4 — Complete benchmark suite
- Run reproducible 5/10/20-qubit circuit grids, warmups, repetitions, and CSV/JSON reports.
- Generate plots and benchmark summary. Never fabricate measurements.

### M5 — Optional extensions (only after M1–M4 work)
- Add small unmerged DFS control.
- Add additional explicitly heuristic priority mode(s).
- Add error-tolerance stopping using the valid frontier bound.
- Explore bounded-memory hybrid DFS without silently dropping contributions.

## 12. Suggested project structure

```text
experiments/best_first_pauli/
├── Project.toml
├── Manifest.toml
├── README.md
├── src/
│   ├── BestFirstPauli.jl
│   ├── encoding.jl
│   ├── gates.jl
│   ├── frontier.jl
│   ├── search.jl
│   ├── circuits.jl
│   ├── baselines.jl
│   └── telemetry.jl
├── scripts/
│   ├── setup.jl
│   ├── run_tests.jl
│   ├── run_benchmarks.jl
│   └── plot_results.jl
├── test/
│   └── runtests.jl
├── results/
│   ├── raw/
│   ├── figures/
│   └── report.md
└── benchmark_config.toml
```

The structure is a suggestion; adapt to the existing repository's conventions.

## 13. Acceptance criteria: definition of done

- [ ] All dependencies installed automatically via documented setup; reproducible environment committed.
- [ ] `BestFirstPauli` runs on supported circuits and produces correct fixed-angle expectation values.
- [ ] Priority-based max-heap expansion, deterministic tie-breaking, and compact path representation implemented.
- [ ] Exact-mode results match official upstream baselines and dense checks on tractable cases.
- [ ] Resource limits are respected; incomplete cases are labeled with a valid frontier error bound.
- [ ] Official upstream `PauliSum` and `VectorPauliSum` benchmarked where supported, with version recorded.
- [ ] Configurations for **5, 10, and 20 qubits** are executed or documented as resource-limited with reproducible evidence.
- [ ] Tests pass, raw data saved, plots generated, and concise report written.
- [ ] No claims of best-first superiority without completed equivalent-work comparisons.
- [ ] All code and commands work on the Intel Core i7-13700HX CPU with 16 GB RAM; no GPU required. The RTX 4060 is optional for later isolated experiments.

## 14. Codex execution instructions

Proceed autonomously through M1–M4. Inspect the actual upstream APIs before coding; adapt constructors and gate implementations to the resolved version. Prefer a small correct implementation over broad but untested support. Run tests after each milestone; fix errors before moving on. Use measured, conservative caps for initial 20-qubit experiments to avoid exhausting the target 16 GB laptop. If long experiments cannot complete within resource limits, produce the partial benchmark data with explicit status rather than pretending success. Finish by summarizing files created, tests run, benchmark commands, results observed, and open limitations.

## References

- Official upstream repository and README: https://github.com/SparqleSim/PauliPropagation.jl
- Upstream releases/version notes: https://github.com/SparqleSim/PauliPropagation.jl/releases
- Pauli propagation framework: https://arxiv.org/abs/2505.21606
