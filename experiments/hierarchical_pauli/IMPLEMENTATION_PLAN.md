# Implementation plan

This isolated experiment reuses the repository's `BestFirstPauli` PauliWord and
gate semantics, whose upstream adapter is verified against PauliPropagation.jl
0.9.2. Maps are sparse dictionaries of local Pauli columns; no dense global
transfer matrix is constructed.

1. Experiment 1: exact sequential columns, map composition, balanced/left-deep
   agreement, and compression telemetry.
2. Experiment 2: the same merged gate kernel with truncation only at block
   boundaries, compared with the gate-wise baseline.
3. Experiment 3: recursive exact maps, a support/column safety cap, and
   hierarchical truncation. All reported end-to-end timings include map build.

The research result is deliberately measurement-driven: a negative compression
or speed result is retained in the reports rather than converted into a claim.
