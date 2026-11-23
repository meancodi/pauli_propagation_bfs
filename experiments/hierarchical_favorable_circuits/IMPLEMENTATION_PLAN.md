# Favorable-circuit study plan

This experiment preserves all existing hierarchical results and uses a new
directory and result namespace. It reuses `BestFirstPauli` gate semantics,
`HierarchicalPauli` exact sparse local maps, and `HierarchicalPauliV2` for
bounded recursive plans and internal-node truncation.

1. Implement deterministic generators for echo, repeated local blocks,
   Clifford sandwiches, perturbed Ising Trotter blocks, and finite-order or
   approximately periodic motifs.
2. Validate every generator against dense references on small widths and check
   Heisenberg ordering, inverse/periodic controls, and exact coefficient maps.
3. Implement a shared paired harness for upstream exact propagation, ordinary
   gate-wise and block-boundary truncation, cached local maps, recursive v2
   maps, repeated-map powers where safe, and explicit on-demand columns.
4. Run a bounded structural screen first, then a moderate campaign with fixed
   seeds, warmups, repeated measurements, durable JSONL records, and resource
   admission checks. Results are written only below
   `results/favorable_v1/`.
5. Generate grouped cross-family plots and a report that ranks positive and
   negative outcomes. A win is accepted only if it survives comparison with
   block-boundary truncation, cached maps, and the optimized upstream backend.

The Windows runner records that memory enforcement is conservative admission
control rather than a hard RSS kill. Unsupported or resource-limited cases are
recorded, not silently replaced.
