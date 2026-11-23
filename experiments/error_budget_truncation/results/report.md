# Error-budgeted coefficient truncation: preliminary research report

## Decision

The bounded representative sweep does not support a claim that the online cost-aware policy improves performance. On the measured 5-, 10-, and 20-qubit cases it is generally comparable to or slower than uniform allocation, while its selection overhead is nonzero. The implementation and benchmark harness are ready for larger repeated studies, but the present evidence is a negative result for adopting a more complicated heuristic.

## What was measured

The run used Julia 1.13.4, PauliPropagation.jl 0.9.2, one Julia thread, seed 11, Ising/random/cancellation-heavy circuits, depths 1 and 2 for 5 and 10 qubits, and depth 1 for a representative 20-qubit run. Budgets were 0, 0.1, 0.01, and 0.001; each measured point had two post-warmup repeats. The raw file contains 600 records and no incomplete records in this run.

At budget 0.1, uniform was the fastest budget policy in the representative records for Ising at all three widths (roughly 0.024 ms, 0.030 ms, and 0.041 ms for 5/10/20 qubits, using the first two-repeat samples as a compact indication). Cost-aware was not consistently better: its current-term-count signal does not create a reproducible runtime advantage over the simple schedule. The cancellation family shows that schedules can consume different amounts of the same allowance, but lower consumption is not itself a runtime win.

The absolute times are small and are useful mainly for relative comparisons; the benchmark records allocated bytes and GC time, but process RSS is not used as a per-case peak-memory certificate. Larger circuits and repeated seeds are needed before making a strong generalization claim.

## Correctness and accounting

The test suite passes 35 assertions. It covers zero-budget exact propagation, cancellation before truncation, deterministic ties, exact ledger conservation, global budget enforcement, dense state-vector agreement, all policies with truncation disabled, fixed-threshold instrumentation, and incomplete resource status.

Every certified checkpoint merges equal Pauli words first. The only discard path adds `sum(abs(c))` to the ledger once. Unused allowance carries forward. Fixed-threshold results are explicitly not labeled with a global certified budget; they are a runtime/error trade-off curve.

## Interpretation

Uniform, front-loaded, and back-loaded schedules are the required simple baselines. The online strategy currently uses only retained-term count and recent growth. In this implementation, selecting terms by a full deterministic sort is measurable overhead, so any adaptive policy must create a sufficiently large propagation saving to pay for selection. The preliminary results do not show that payoff.

Certified bounds are conservative: in several records the observed expectation error is zero while the cumulative discarded mass is positive. This is expected because the coefficient-mass bound is state-independent and does not exploit cancellations in the final observable.

## Reproduction and artifacts

See [`README.md`](../README.md), [`IMPLEMENTATION_PLAN.md`](../IMPLEMENTATION_PLAN.md), and the pinned [`Manifest.toml`](../Manifest.toml). Raw JSONL data is in [`results/raw/benchmark.jsonl`](raw/benchmark.jsonl), with summary tables in [`summary.md`](summary.md) and figures in [`figures/`](figures/).

The full configured sweep (depths 1, 2, and 4; seeds 11 and 29) should be run before publication-quality claims. Do not replace the bounded result above with fabricated results if that sweep exceeds the available memory or runtime; retain incomplete records and report them as such.
