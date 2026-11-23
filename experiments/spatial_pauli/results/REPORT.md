# Spatially partitioned Pauli propagation: exact experiment

## Research objective

The experiment tests whether qubit-space partitioning and lazy tensor-product factors improve Pauli propagation. It is independent of the earlier gate-wise hierarchical transfer-map experiment. The result is mixed: one narrow cross-coupled circuit benefits, while most measured cases are substantially slower.

## Method and correctness

S0 uses a conventional global sparse Pauli dictionary. S1 uses two contiguous qubit groups. S2 uses contiguous two-qubit leaves in a balanced spatial tree. A component stores a scalar times one sparse local Pauli sum per leaf. Thus a product of two two-term factors uses two local sums, without expanding four global terms. Every circuit uses the same gates, parameters, observable, and computational zero state. Gates are propagated in reverse stored circuit order, as required for Heisenberg conjugation.

Local gates touch only their leaf factor. Cross-leaf Pauli rotations split touched factors by commutation parity; the odd-parity sectors receive the exact sine branch, with Pauli product phases preserved. Cross-leaf CNOTs split the control and target factors by Pauli letter and apply the generic tested Clifford conjugation kernel. The tree identifies the lowest common ancestor of affected leaves. Exact aggregation normalizes and combines identical products and sums final factors when all earlier factors match. No global operator expansion occurs in the propagation path. A full expansion exists only for testing. No approximate compression was enabled because the exact method did not show broad benefit.

The tests compare every coefficient against S0 for small deterministic and randomized circuits, including non-power-of-two sizes, identities, cancellation, local and cross-leaf rotations, CNOTs, and noncommuting gate sequences. Three 4-qubit cases were also checked against the upstream `PauliPropagation.jl` reference and independent dense-state expectations. **157/157 tests passed.** The largest observed benchmark zero-state expectation difference against S0 was **4.44×10⁻¹⁶**.

## Benchmark method and scope

The corrected benchmark records, identified by the `CPU` column, use Julia 1.13.4 on Windows, CPU name `raptorlake`, and the thread count recorded in each CSV. Fixed seeds generate identical circuits for all methods. Each method gets one warm-up, three timed repetitions, and a separate memory-profile pass. Reported time is the median end-to-end time including construction and aggregation; min and max are retained. Allocated bytes are measured in the timed runs. Peak live bytes are the maximum Julia `summarysize` of the live dictionary or spatial operator during the separate profile; they are **not** operating-system RSS. The temporary contribution count is recorded separately. A 5,000-component, 100,000-local-term, and 20,000-contribution limit stops work rather than falling back.

The corrected campaign covers all five families, four observables, and depths 1 and 2 at 4 and 8 qubits; all five families at 12 qubits and depth 1; pair and Clifford controls at 16 qubits and depth 1; and a separate fixed-depth 8-qubit weak-coupling sweep with 0–3 boundary interactions. This is **333 distinct method records** and **222 paired S1/S2 comparisons against S0**. Earlier CSVs lacking the `CPU` field include an exploratory campaign whose timing was contaminated by per-gate memory inspection. They are preserved as raw history and excluded from the figures and conclusions. One exploratory 12-qubit cross-coupled configuration was stopped after excessive runtime; it is not counted as a successful or failed corrected result. The corrected grid is therefore bounded rather than exhaustive at 12–16 qubits.

## Results

Only **4 of 222** paired spatial runs were faster than S0. The median spatial/S0 runtime ratio was **19.6×**; the median estimated live-memory ratio was **2.80×**. Estimated memory was lower in **57 of 222** pairs. These counts treat individual configurations as observations, not independent statistical replicates.

The clearest advantage was the 8-qubit, depth-2 random cross-partition circuit with a central ZZ observable. A repeated corrected run measured S0 at **7.685 ms**, S1 at **2.222 ms** (71% faster), and S2 at **4.078 ms** (47% faster). Estimated live representation storage was **409,712 B** for S0, **13,992 B** for S1, and **17,864 B** for S2. Another corrected run of the same seeded case also found both spatial methods faster, though the relative order of S1 and S2 varied. This supports a **narrow, reproducible advantage**, not a general speedup.

The worst measured runtime ratio was S2 on the 12-qubit independent-pair circuit, depth 1, with a local-sum observable: **11.812 ms versus 0.0751 ms** for S0, or **157× slower**. Its estimated live storage was **30.2×** S0. The smallest memory ratio was S2 on an 8-qubit cross-partition circuit, depth 1, with a multi-block observable: **2.58%** of S0's estimated live storage, while runtime was **1.36× slower**. At 16 qubits, the sampled pair and Clifford cases all favored S0 in runtime. S2 did not consistently improve on S1.

Figures from corrected runs: [total runtime](figures/runtime.png), [estimated live memory](figures/memory.png), [factorized component growth](figures/rank_growth.png), [coupling sensitivity](figures/coupling_sensitivity.png), and [percentage improvement against S0](figures/improvement.png). The improvement figure has a zero baseline; positive values denote improvements and negative values regressions. No accuracy-versus-runtime figure is provided because no approximate method was run. Raw records are in [`raw/`](raw/).

## Interpretation and limitations

The best case shows that parity-sector updates can avoid a costly global sparse expansion. But factorized components, normalization, dictionary allocations, and tree traversal often outweigh that benefit. The independent-pair and Clifford controls are particularly unfavorable because S0 remains compact and the spatial bookkeeping has little to save. Cross-partition gates can rapidly increase the number and size of local factors; exact aggregation is not an optimal tensor-rank reduction. The tree's flat leaf-factor representation also incurs overhead and does not exploit every possible subtree reuse. The coupling sweep at one depth and observable did not establish a monotonic trend; its observable and gates sometimes commute, so it should not be used to infer a general coupling law.

The memory figure measures retained Julia objects, not transient allocator peaks or RSS. The diagnostic `peak_contributions` counts intermediate components but is not a byte estimate. The benchmark does not enforce a preemptive per-case wall-clock timeout, and it does not report a certified error bound because the implementation is exact. The 12-qubit depth-2 and broader 16-qubit grids were not completed under the corrected timing method. The raw files and benchmark script retain the exact configurations and limits for extension.

## Conclusion

Spatial partitioning can reduce both runtime and estimated representation memory for a specific cross-coupled case, but **does not provide a general improvement** in this campaign. Conventional sparse propagation is superior for most local, independent-pair, and Clifford-dominated controls measured here. Two partitions were usually cheaper than the multilevel tree, so the tree is not justified as the default. Approximate compression remains a research option only if a future workload demonstrates sustained low rank and a useful accuracy–cost tradeoff; the present exact results do not justify implementing it now.
