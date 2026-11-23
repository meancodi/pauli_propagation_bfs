# Spatial Pauli V2: exact qubit-wise hierarchy

## Objective and V1 audit

V2 asks whether a genuine recursive qubit-space operator representation improves on the corrected flat factorization and conventional sparse propagation. It does not extend the gate-wise `HierarchicalPauli` experiments.

V1 had four material problems. Reversed cross-leaf CNOTs used spatial block order as control/target order; a new regression initially failed 48 coefficient checks. Recursive interval splitting produced misaligned leaves, including singleton leaves for a requested pair partition. Exact aggregation only combined the last factor, and normalization could pivot on a tiny coefficient. Finally, the initial benchmark included `summarysize` in the timed gate loop and used incomplete plotting keys. V1 reports, CSVs, and figures remain untouched; V2 uses separate raw records and figures.

## Exact methods

S0 is a global sparse Pauli dictionary. `S0_optimized` uses the installed PauliPropagation.jl `PauliSum` backend with zero coefficient threshold, unlimited Pauli weight, and one thread. S1 is an aligned two-block flat factorization. S2 is an aligned flat vector of two-qubit leaf factors. S3 is a recursive operator DAG: each leaf stores a sparse local sum, while each internal node stores an exact sum or tensor product of child expressions. Local gates reuse unaffected sibling subtrees. Boundary rotations split internal expressions by commutation parity; boundary CNOTs project control and target letters and recombine at the owning ancestor. This is actual internal-node algebra, rather than a flat vector with a decorative tree.

The representation supports exact addition, scalar multiplication, and tensor construction. It preserves stored gate order in reverse for Heisenberg propagation. The full global expansion routine is used only for small-system tests. The recursive implementation may use complex intermediate coefficients for Pauli product phases; final Hermitian coefficients and expectations agree with the real reference to Float64 tolerance. No approximate compression or S4 switching method was implemented.

## Correctness

The final suite passed **2,416/2,416** checks: CNOTs in both directions and all I/X/Y/Z control-target letters, dense CNOT matrix conjugation, aligned partitions for 4–16 qubits and leaf sizes 1–4, 64-bit masks, X/Y/Z and multi-site rotation phases, non-power-of-two S3 trees, periodic endpoint gates, randomized full-operator comparisons, expectation values, cancellation, and invalid-input rejection. Details are in [CORRECTNESS.md](CORRECTNESS.md). The largest benchmark expectation error against S0 was **1.11×10⁻¹⁵**. Small-system tests compare every Pauli coefficient, not just expectations.

The small-system Schmidt-rank diagnostic forms a coefficient matrix only after exact propagation. Across the 4|4 cut, the sampled 8-qubit pair circuit had rank 1, a depth-1 cross-partition factorized case had rank 2, and a depth-2 cross-partition ZZ case had rank 8. S2's peak component counts matched these three ranks in this sample, but component count is **not** generally mathematical rank. See `raw/v2_rank_20261122_230047.csv`.

## Measurement and coverage

[BENCHMARK_METHOD.md](BENCHMARK_METHOD.md) defines seeds, circuit fingerprints, all configuration fields, timing categories, resource limits, and memory units. Each successful method has a warm-up, five timed repetitions, and a separate memory-profile pass. The V2 plots use the latest pilot, partially completed targeted, completed 10-seed, final, and coupling files. They contain **373 distinct completed method records** and **75 configurations with an S0 record**; 74 have matched S2 and S3 records. A stopped exploratory targeted campaign preserved 78 rows, leaving one configuration missing S2/S3 and 64 planned configurations unattempted. An empty final CSV from a corrected script error is preserved but excluded. No failed or incomplete case was counted as a spatial success. A separate four-case, 20-row diagnostic file records per-node expression arity.

`peak_retained_bytes` is Julia `Base.summarysize` of the live representation, not RSS or peak temporary allocation. Allocated bytes are the maximum across five timed runs. The `peak_intermediate` field counts expressions or components, not bytes. The campaign did not enforce a preemptive per-case timeout; the interrupted targeted sweep is reported explicitly. Initialization and propagation are separated for S0–S3. For the upstream backend, its combined API call is timed end to end, and its subphases cannot be separated fairly. S1/S2 gate time includes aggregation, with a separate measured aggregation value and a non-overlapping gate-excluding-aggregation value. S3 internal aggregation events are counted, but their time is included in its gate phase.

## Optimization result

Aggregation dominated approximately 95% of S1/S2 gate time in the first V2 pilot. A single-component fast path skips signature building and multi-component aggregation while retaining exact zero removal and resource checks. On 12 matched pilot configurations, median S1 runtime fell to **0.346×** and S2 to **0.315×** their pre-change values; allocated bytes fell to **0.225×** and **0.223×**. All tests remained correct. This is a within-V2 before/after comparison using identical fingerprints, documented in [OPTIMIZATIONS.md](OPTIMIZATIONS.md). The historical V1 median slowdown cannot be compared numerically with V2 because circuits, seeds, partitions, and timing semantics changed.

## Performance results

Across the 75 matched configurations, median total-time ratios against S0 were **1.06×** for `S0_optimized`, **5.64×** for S1, **7.93×** for S2, and **3.99×** for S3. S3 beat S0 in **12 of 74** comparisons. The median S3/S2 ratio was **0.538×**, and S3 beat S2 in **58 of 74** comparisons. Thus the recursive structure delivers a reproducible gain over flat multi-leaf factorization, but conventional sparse propagation remains faster in most measured cases. S3 beat the upstream backend in only **6 of 74** matched cases.

Median retained-size ratios against S0 were **0.451×** for S1, **0.741×** for S2, and **0.943×** for S3. Median allocated-byte ratios were **1.33×**, **1.99×**, and **1.98×**, respectively. The recursive DAG tends to save retained storage in some cases while allocating more during updates. S3's median retained size was **1.20×** S2's, despite its runtime advantage over S2.

The best S3/S0 time result was the 8-qubit, depth-2 random cross-partition ZZ case: **2.436 ms versus 25.019 ms**, a **90.3% improvement**, with retained size **28,744 B versus 1,638,512 B**. A second seed of the same structure also favored S3 (**1.756 ms versus 6.701 ms**). The upstream backend took **2.655 ms** and **1.553 ms** on those seeds, so S3 was competitive with it only in the first. These two seeds are evidence for a narrow opportunity, not a broad scaling claim.

The worst measured S2/S0 runtime ratio was a periodic 8-qubit Ising ZZ control: **1.335 ms versus 0.0233 ms**, or **57.3× slower**. The 10-seed study of factorized multi-term observables found S3 faster than S2 on all 20 seeded cases, with median ratios **0.537×** for independent pairs and **0.222×** for cross-partition circuits. S3 still had median slowdowns against S0 of **1.16×** and **4.15×** in those two families. The fixed-depth weak-coupling sweep at 0–3 boundary interactions did not show a reliable monotonic performance trend.

Figures: [runtime](figures/runtime.png), [retained size](figures/memory.png), [allocations](figures/allocations.png), [percentage improvement](figures/improvement.png), [component growth](figures/component_growth.png), [coupling sensitivity](figures/coupling_sensitivity.png), [S2 versus S3](figures/s2_vs_s3.png), and [construction versus propagation](figures/construction_vs_propagation.png). The improvement plot has a zero baseline, positive for improvement and negative for regression. It shows one measured seed per family and method for readability; paired statistics above use all valid records.

## Adaptive feasibility and limitations

A cheap selector based only on initial observable term count and boundary-gate frequency is not supported by these measurements. All 20 seeded factorized-observable circuits began with four terms; S3 beat S0 on only four of them, and the cross family was usually slower even though two deeper cross-coupled ZZ cases were much faster. Selecting the latter from completed outcomes would be oracle tuning. Conversion costs and potential representation oscillation have not been measured, so an S4 switcher is not justified.

S3 retains the unchanged subtree and performs useful internal exact aggregation, but `rsum` can perform many structural comparisons and allocations. Its internal sum arity is not optimal tensor rank. Cross-child gates may grow expression size rapidly. Per-node arities and physical gate crossing levels are available in the diagnostic file, while the main campaign records reuse and internal-aggregation counts. Process RSS, peak temporary bytes, exact S3 aggregation time, and enforced wall-clock case limits remain open measurement gaps. The broad 12–16-qubit and 10-seed depth-2 cross-family grids were intentionally bounded after the costly targeted run.

## Conclusion

The V2 repair makes the flat methods correct for reversed CNOTs and aligned partitions. The single-component optimization measurably reduces S1/S2 overhead. S3 is a genuine recursive representation and is **consistently better than S2** in the 10-seed study, but it is **not a generally better replacement for S0 or the upstream backend**. Further S3 work should focus on bounded internal aggregation and structural sharing for the specific cases where the global sparse map becomes large. Approximate compression and adaptive switching should wait for stronger, prospectively predictable opportunities.
