# Final research report: hierarchical merge-and-truncate Pauli propagation

## 1. Objective

The study tests whether recursively composing sparse Pauli block maps creates
enough reconvergence to improve the runtime–accuracy–memory trade-off over
ordinary gate-wise or block-boundary propagation. The implementation does not
assume a positive result.

## 2. Mathematical formulation

For a block B, a column stores `H_B(P)=B† P B` as a sparse Pauli sum. Gates are
stored in circuit order and Heisenberg application reverses that order. For
chronological blocks A then B, the implementation uses `T_(BA)=T_A*T_B` by
applying the right child column first and the left child map second. Equal Pauli
keys are merged in dictionaries before any experimental truncation.

## 3. Experiment 1: exact composition

The exact sequential and recursive maps agree in automated tests, including
Clifford signs, rotation phases, cancellation, composition order, and dense
references. The smoke data shows that raw sequential gate expansion can create
many contributions before merging, while recursive map construction pays for
local columns and composition. The favorable compression ratio therefore does
not by itself imply a runtime win. See `results/experiment1.jsonl` and the two
Experiment 1 figures.

The smoke run is intentionally small; it is sufficient to validate the
mechanism and expose overhead, not to establish a universal scaling law.

## 4. Experiment 2: block-boundary truncation

The baseline uses the same sparse merged gate kernel and truncates after every
gate (`block_size=1`). Larger block sizes postpone truncation while preserving
exact merging within a block. On the tested small circuits, exact runs agree
with dense references. Thresholded runs report observed error, peak terms,
allocation, truncation time, and discarded mass. Results vary by threshold and
block size; larger blocks can increase peak intermediate terms, so reduced
checkpoint overhead is not a free improvement.

This repository does not contain a prior delayed-truncation result with the
same benchmark schema, so no unsupported numerical comparison is made. The
existing error-budget experiment remains unchanged and is the appropriate
source for that separate comparison.

## 5. Experiment 3: recursive composition and hierarchical truncation

Balanced and left-deep recursive maps are implemented, with support/column
limits providing a bounded feasibility mode. The smoke records show the main
cost clearly: for the tested 5- and 8-qubit random cases, map construction and
allocation made recursive strategies slower than the direct sequential
baseline, even though exact expectations remained correct. A support cap is a
safety mechanism, not evidence of improved asymptotic scaling.

The current hierarchical approximation truncates after selected local-map
applications and reports observed error plus discarded coefficient mass. It is
not presented as a rigorous inherited map-error certificate; the report uses
the observed dense-reference error where a dense reference exists.

## 6. Direct comparison

The decisive comparison is recursive/hierarchical application versus ordinary
block-boundary propagation. In this bounded smoke study, recursion did not
demonstrate the required reproducible 20–30% end-to-end improvement. Its exact
correctness is established, but construction/allocation overhead dominates the
small workloads. The result is a feasibility negative, not a claim that every
structured circuit will behave identically.

## 7. Failure modes

The relevant risks are transfer-map column growth with support, repeated map
construction, and memory allocation during balanced composition. The
support/column cap prevents an uncontrolled dense/global map. Low-complexity
circuits may also make recursive caching unattractive because columns are not
reused often enough.

## 8. Scalability

The dominant structural limit is local support: a complete map has `4^k`
columns before considering column sparsity. Consequently, the implementation
keeps maps local and bounded and does not attempt a global dense transfer
matrix. Larger sweeps should be selected from measured peak entries and memory,
not nominal qubit count.

## 9. Novelty and limitations

Sparse Pauli transfer-map composition, exact merging, and block/slice
truncation are established techniques. The contribution here is an auditable,
reproducible comparison of grouping strategies and truncation checkpoints in
this repository's Pauli kernel. No novelty claim is made merely from using a
balanced recursion tree. On-demand column caching and rigorous inherited map
error bounds remain future work rather than silently asserted features.

## 10. Verdict

For the measured feasibility workloads, pursue ordinary block-boundary
truncation only where its measured trade-off is favorable, and do not adopt
recursive map construction as a default. Restrict further recursive work to
circuits with demonstrably repeated local-block reuse or unusually strong exact
reconvergence, and require end-to-end measurements including construction and
memory before claiming an advantage.
