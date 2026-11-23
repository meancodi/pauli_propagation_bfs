# Spatial Pauli experiment

V2 continues this experiment with corrected CNOT orientation and partitions, an exact recursive S3 operator, and separate results in [`results/v2/REPORT.md`](results/v2/REPORT.md). Run the V2 benchmark with `julia --project=experiments/hierarchical_pauli experiments/spatial_pauli/benchmarks/run_v2.jl` from the repository root.

This experiment factors an observable into a sum of tensor products of sparse local Pauli sums. It is independent of the gate-wise hierarchical experiment. Gates are applied in reverse circuit order for Heisenberg propagation. The reference is a conventional global sparse dictionary using the same tested Pauli-word conjugation kernel.

## Exact boundary update

For a cross-block Pauli rotation, each touched local factor is split into its commuting and anticommuting parts relative to the local piece of the generator. The Cartesian product of **parity sectors** is formed, never the Cartesian product of all Pauli terms. Odd global parity sectors receive the cosine branch and a generator-multiplied sine branch. Local Pauli product phases are combined to give the correct real sign. A cross-block CNOT partitions the two touched factors by their Pauli letters at the control and target, at most 16 sector pairs. The generic Clifford conjugation kernel supplies the mapped letters and sign. Gates inside a leaf update that leaf's sparse sum only. The tree defines contiguous leaves and spatial ancestry; a flat vector of leaf factors is a factored representation at every internal cut. Tree nodes do not store transfer matrices or copies of operator data.

Exact aggregation normalizes factors and combines identical products. It also sums final factors when all earlier factors match. This is valid rank reduction, but it is not an optimal low-rank factorization. No approximate compression is implemented.

Run tests with `julia --project=../hierarchical_pauli test/runtests.jl` from this directory. The existing project supplies the installed generic kernels; no hierarchical code is imported. Benchmarks can be run with `julia --project=../hierarchical_pauli benchmarks/run_benchmarks.jl`. Results are appended to a timestamped CSV in `results/raw/`.

Limits on component count, local terms, and intermediate contributions are enforced after each update or during branch construction. They stop a case with a recorded status; there is no fallback. The benchmark grid is deliberately bounded. Peak live representation is a Julia `summarysize` estimate, not process RSS. Wall-clock duration is recorded but not preemptively enforced.
