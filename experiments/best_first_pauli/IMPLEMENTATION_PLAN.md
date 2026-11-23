# Implementation decisions

1. Resolve official registered PauliPropagation and inspect its actual source/API.
2. Implement an independent UInt64 symplectic Pauli representation, real signed
   rotation/Clifford kernel, deterministic binary max-heap and bounded search.
3. Test dense conjugations and state evolution, upstream agreement, ordering,
   interruption bounds and reproducibility before benchmarks.
4. Run both circuit families at 5/10/20 qubits, multiple layers, seeds and two
   observables. Use isolated worker processes with external RSS/runtime guards.
5. Save repeat-level CSV, configuration/metadata, SVG figures and a Markdown report.

The repository initially had only a README and the supplied design. This experiment
is isolated under experiments/best_first_pauli. Gates are stored in Schrodinger
order and indexed from length down to 1 (the design pseudocode ambiguously combines
reversed storage with descending indexing). Fixed-angle Float64 Hermitian Paulis
have real coefficients; complex multiplication phases are resolved into real signs.
No path merging, coefficient truncation, GPU or hybrid DFS is introduced.

## DFS and equal-budget extension

1. Add an explicit LIFO backend using the same `Path`, `expand`, terminal overlap,
   compensated sums, limits, and frontier-bound invariant as best-first.
2. Define deterministic DFS order by pushing sine second in traversal order: the
   implementation pushes children in reverse, so cosine/unchanged is visited first.
3. Validate exact equality, partial bounds, frontier retention, deterministic order,
   and genuinely equal expansion counts before collecting measurements.
4. Run node-budget, runtime-budget, and certified-tolerance experiments in isolated
   single-thread workers. Compute exact references separately with the upstream
   merging backend and label cases without a completed reference as unavailable.
5. Store all extension outputs under `results_dfs_vs_bestfirst`, preserving `results`.
