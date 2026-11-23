# Local Coefficient Blocks: implementation plan

## Scope

This isolated experiment implements a fixed-region Local Coefficient Block (LCB)
representation for exact Pauli propagation, followed by regrouping and a small
adaptive policy only if the fixed-region measurements justify it. Existing
experiments and the upstream PauliPropagation.jl package are not modified.

## Repository findings

- `BestFirstPauli.PauliWord` uses a symplectic `UInt64` pair `(x,z)`, with qubit 1
  in the low bit; `Y` has both bits set. Coefficients are `Float64` and the
  expansion kernels resolve Pauli-product phases into real signs.
- `BestFirstPauli.expand` is the independent local gate kernel used by the
  existing tests and by this experiment for parity with the established
  Heisenberg convention.
- The official baseline is `PauliPropagation` 0.9.2 through its public
  `PauliSum`, `VectorPauliSum`, and `propagate` APIs. Its adapter is already
  validated in `experiments/best_first_pauli/src/baselines.jl`.
- `experiments/hierarchical_pauli` groups transfer-map columns and composes
  local maps. That is related partitioning work, but it does not store an
  evolving Pauli sum as outside-key -> local coefficient block, so it is not an
  equivalent implementation.

## Design

For active qubits `S`, the outside key is `(x & ~mask, z & ~mask)`. A local
configuration is packed in base four, with digits `I,X,Z,Y` in active-qubit
order. Dense blocks contain `4^k` values; sparse blocks contain sorted local
indices and values. Gate kernels operate directly on local indices when a gate
support is wholly inside `S`. Gates crossing the boundary use a safe conversion
to ordinary terms followed by regrouping.

The implementation supports conversion from/to `Dict{PauliWord,Float64}` and
convenience conversion from upstream sums. It exposes regrouping, iteration,
zero removal, exact Clifford propagation, and exact one-/two-qubit Pauli
rotations.

## Execution order

1. Implement and test dense/sparse fixed-region LCB.
2. Run small correctness checks against the independent expansion kernel and
   upstream expectation values.
3. Run reproducible synthetic locality/occupancy benchmarks.
4. Add regrouping and measure switching break-even points.
5. Add a conservative look-ahead adaptive wrapper only if fixed-region results
   show a favorable regime.

## Resource policy

Primary runs use one Julia thread, Float64, fixed seeds, warm-up before timing,
JSONL output, and a conservative term cap. Million-term cases are attempted
only when the measured host memory allows them; skipped/OOM/timeout cases are
recorded rather than retried indefinitely.
