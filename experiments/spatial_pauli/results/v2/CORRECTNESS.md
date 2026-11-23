# V2 correctness

The V1 flat representation failed for CNOTs whose control was in a later leaf than the target. The V2 regression first reproduced 48 failing coefficient checks. `cross_cnot!` now finds control and target leaves independently and keeps gate orientation throughout letter projection and reconstruction.

The partition constructor now builds contiguous fixed-size leaves first and then balances the leaf sequence. For 12 qubits and leaf size 2, the leaves are exactly `(1,2), (3,4), (5,6), (7,8), (9,10), (11,12)`. The 64-bit full-support mask has an explicit case. Gate support is validated before propagation.

Flat exact aggregation combines common factors at any leaf position, for at most the configured number of passes. Normalization chooses the largest-magnitude nonzero coefficient, breaking ties by Pauli bit pattern. It rejects nonfinite input and never drops a small nonzero coefficient by threshold. Detected underflow in normalization or local scaling and detected normalization overflow raise explicit errors. Other extreme floating-point products can still underflow or overflow; exact mode provides no arbitrary-precision guarantee.

S3 stores a recursive expression DAG: leaves hold local sparse sums; internal nodes hold tensor products or sums of expressions over child subtrees. Local gates reuse the unchanged sibling expression. Cross-child rotations split expressions into commuting and anticommuting sectors, multiply the generator with complex Pauli phases, and recombine exactly. Cross-child CNOTs project the control and target letters, transform them with the tested Clifford kernel, and recombine at their ancestor. The global expansion routine is used only in tests.

Run `julia --project=experiments/hierarchical_pauli experiments/spatial_pauli/test/runtests.jl` from the repository root. The suite checks full coefficients against S0 for randomized 4, 5, 6, and 8-qubit circuits, both CNOT orientations, X/Y/Z and multi-site rotations, periodic endpoints, non-power-of-two trees, identity factors, cancellation, and dense matrix conjugation. Float64 coefficient tolerance is `1e-10` for circuits and `1e-12` for single-gate checks.
