# V2 progress

- Reproduced reversed CNOT bug with 48 failing coefficient checks, then fixed orientation.
- Aligned fixed-size leaves, validated masks and gate support, stabilized normalization, and extended exact aggregation across factor positions.
- Implemented S3 as a recursive operator expression DAG with internal sums, tensor products, and subtree reuse.
- Full coefficient and dense checks pass; pilot has 60 completed method records.
- Benchmarked the single-component S1/S2 fast path against a matched pre-change pilot; see `OPTIMIZATIONS.md`.
- Final suite: 2,416/2,416 checks passed. Completed 10-seed pair and cross-partition studies, a bounded 8/12/16-qubit campaign, a coupling sweep, and a four-case Schmidt-rank diagnostic.
- Generated eight V2 figures and `REPORT.md`. S3 beat S2 in 58/74 matched comparisons but beat S0 in 12/74. A costly targeted sweep was stopped after 78 incrementally saved rows; the report counts its incomplete coverage.
