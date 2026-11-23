# V2 optimization log

## Exact aggregation across every position

- Affected: `aggregate!` in `src/SpatialPauli.jl`.
- Correctness: complete V2 tests, including cancellation and full coefficient comparisons.
- Effect: enables `A₁⊗B + A₂⊗B = (A₁+A₂)⊗B` for any leaf position. This is an exact simplification, not a minimum-rank search. The bounded passes can add work when no common factors exist.

## Single-component fast path

- Affected: `finalize_components!` after gate propagation.
- Correctness: full V2 suite passed after the change; maximum pilot expectation error was `2.22e-16`.
- Comparison: 12 identical pilot configurations in `raw/v2_pilot_20261122_224126.csv` (before) and `raw/v2_pilot_20261122_224348.csv` (after), matched by complete circuit fingerprint and method.
- S1 median runtime after/before: **0.346**; allocated bytes after/before: **0.225**.
- S2 median runtime after/before: **0.315**; allocated bytes after/before: **0.223**.
- S3 median runtime after/before: **0.983**; allocated bytes ratio: **1.000**, as expected because S3 does not use this path.

When there is one factorized component, there are no other components to merge. The fast path removes exact zeros, checks finite coefficients and limits, and skips normalization and signature construction. If more components appear later, full exact aggregation resumes. This change was retained because it improved both runtime and allocation on matched pilot cases.
