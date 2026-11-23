# Implementation plan and API decisions

1. Reuse the existing independent `BestFirstPauli` PauliWord encoding, circuit generators, dense references, and upstream adapter. The official dependency is PauliPropagation.jl 0.9.2 in the isolated environment.
2. Implement a separate dictionary-based merged Pauli sum and a straightforward stable sort for coefficient selection. Equal magnitudes use `(x,z)` word keys; this makes decisions reproducible and keeps selection overhead measurable.
3. Propagate in the same Heisenberg gate order as the existing upstream comparison, merge before truncation, and make checkpoint allocation the only certified discard path.
4. Implement uniform, linear front-loaded, linear back-loaded, and a small online cost-aware policy. Unused budget is carried forward through the remaining-checkpoint allocation rule.
5. Test exact propagation, cancellation, ledger conservation, deterministic ties, budget limits, independent dense agreement, and incomplete resource status before benchmarking.

The upstream source is not modified. The upstream fixed-threshold implementation remains available through the existing adapter, while the local fixed-threshold baseline additionally records discarded mass for trade-off plots.
