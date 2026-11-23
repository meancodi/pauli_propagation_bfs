# Experiment 1: fixed-region LCB

## Result

The exact fixed-region implementation is correct, but it is not a general
replacement for the upstream backends. Dense blocks show a useful specialized
kernel regime for repeated gates inside the active region; sparse blocks do not
recover the expected advantage at `k=2`.

## Implementation

`DenseLCB` and `SparseLCB` store `outside_key => local coefficients`. The
outside key masks the active qubits from the repository's `(x,z)::UInt64`
encoding. Local indices are base four (`I,X,Z,Y`) and use no strings. Local
Cliffords and one-/two-qubit Pauli rotations call the established exact
`BestFirstPauli.expand` sign kernel directly on local words. Boundary gates use
a safe ordinary-term conversion followed by regrouping.

## Correctness

Tests cover all 16 local indices for noncontiguous active qubits, duplicate
accumulation, cancellation, H/S/CNOT, a two-qubit rotation, crossing gates,
regrouping, and round-trip conversion. Circuit comparisons in
`results/upstream_circuits.jsonl` agree in expectation to approximately
`1e-14` with both upstream backends.

## Measurements

The synthetic campaign uses fixed seeds, Float64, one Julia thread, three
timed repetitions, 1k/10k/50k/100k requested terms, and requested occupancies
5/15/30/60/90%. Actual occupancy is recorded because group rounding changes it.
It uses four repeated local gates and a random/crossing control.

At 50k terms, local-gate medians include the following representative points:

| Requested occupancy | Dictionary | Dense LCB | Sparse LCB |
|---:|---:|---:|---:|
| 5% | 29.9 ms | 41.5 ms | 51.3 ms |
| 15% | 31.2 ms | 15.2 ms | 25.6 ms |
| 90% | 29.9 ms | 33.1 ms | 60.1 ms |

The exact raw records are in `results/synthetic.jsonl`; plots are under
`results/figures/`. The 15% point is a real favorable regime, but it is not
monotonic or universal: dense clearing and local expansion costs dominate at
very low occupancy, while the dictionary wins again in some high-occupancy
cases.

For random/crossing gates, dense and sparse LCB are slower because they scan,
materialize, and regroup the full term set. This is the central limitation of
fixed active regions.

## Upstream comparison

The small full-circuit campaign covers Ising/random, widths 5/8, depths 4/8,
and two seeds. Median elapsed times across these cases were approximately:
`VectorPauliSum` 1.1 ms, `PauliSum` 2.3 ms, dense LCB 8.5 ms, sparse LCB
15.4 ms. These include construction and all gates, so they are the relevant
end-to-end comparison. LCB did not beat either official backend in this
campaign.
