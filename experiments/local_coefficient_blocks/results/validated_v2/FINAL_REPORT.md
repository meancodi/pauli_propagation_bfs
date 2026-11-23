# Corrected LCB benchmark report (v2)

## Executive conclusion

DenseLCB does not establish a reproducible, practically meaningful
end-to-end speedup over optimized `VectorPauliSum` when active-region changes
and representation costs are included. It can be a specialized local kernel:
in the zero-regrouping control, some 10k-term, 60–90% occupancy cases became
faster after sufficiently many local gates. That advantage was not stable in
the alternating-support workload and was absent in the completed 50k-term
case. The recommendation is to discontinue general-purpose/adaptive LCB
development and retain the implementation only as a controlled microbenchmark
or specialized-kernel reference.

## 1. Defects corrected

The original generator sampled unrestricted 64-bit masks. The v2 generator
enumerates outside configurations in base four over only the complement of the
declared active region, under the exact bound `4^(n-k)`. It handles `n=64`
without shifting by 64, validates all emitted masks, records per-block
occupancy, and rejects infeasible requests instead of looping indefinitely.

The previous regrouping experiment changed the active region to `[3,4]` while
continuing to apply gates on `[1,2]`. V2 alternates matching gate supports:
`[1,2]` uses H(1), CNOT(1,2), and an XZ rotation on (1,2); `[3,4]` uses the
equivalent gates on (3,4). The corrected local-only runs report zero crossing
fallbacks.

The old allocation measurement propagated an already-mutated representation.
V2 uses independent copies from the same immutable input for timing,
allocation, and correctness checks. Construction-inclusive and steady-state
timings are separate. JSONL output is append/resume safe.

The old upstream plot connected unrelated widths, depths, families, and seeds.
V2 plots use paired configuration keys and separate mode/width/occupancy
figures.

## 2. Mathematical correctness and upstream path

The LCB continues to use the repository's `(x,z)::UInt64` encoding, low-bit
qubit convention, exact Float64 coefficients, and established Clifford/
rotation sign kernel. V2 tests cover n = 4, 8, 16, 32, and 64, noncontiguous
regions, mask bounds, occupancy feasibility, support matching, cancellation,
and matched DenseLCB/SparseLCB versus `VectorPauliSum` expectations.

The official path is the public PauliPropagation.jl 0.9.2 API:

```julia
PP.propagate(gates, sum, theta;
    min_abs_coeff=0.0, max_weight=Inf, thread=false)
```

The benchmark constructs actual `PP.VectorPauliSum` and `PP.PauliSum` objects,
uses Float64 coefficients, disables truncation, and uses the same reversed
gate-list convention as the existing verified adapter. The primary comparison
is therefore against the optimized upstream propagation path, not the old
custom dictionary control.

## 3. Corrected campaign

The durable file `break_even.jsonl` contains 416 completed backend rows:

| Requested terms | Completed cases | Repetitions |
|---:|---:|---:|
| 1,000 | 48 | 7 |
| 10,000 | 48 | 7 |
| 50,000 | 8 | 5 |

The completed cases cover requested occupancy 15%, 60%, and 90%; primitive
gate counts 1, 4, 16, and 64 per segment; 2 and 8 segments; and both
zero-regrouping and correctly alternating-support modes. A case's actual mean
occupancy is recorded separately; for example, requested 15% gives 12.5%
with `k=2` because local occupancy is quantized in units of 1/16.

For each case, the stored speedup is

\[
S=T_{\mathrm{VectorPauliSum}}/T_{\mathrm{DenseLCB}}.
\]

Thus `S > 1` means DenseLCB is faster. The raw data stores median, MAD,
minimum, and maximum for each timing category.

### Representative results

At 10,000 terms, 8 segments, zero regrouping, and 90% requested occupancy,
DenseLCB's steady-state speedup was about 1.21 at 16 gates/segment and 1.05
at 64 gates/segment. Including construction, the corresponding speedups were
about 1.12 and 1.08. This is the strongest reproducible favorable regime in
the completed matrix.

The same 10,000-term 90% occupancy alternating-support cases gave only about
1.08 at 16 gates/segment and 1.01 at 64 gates/segment end-to-end, while the
steady-state ratios were below one. At 10,000 terms and 60% occupancy, the
alternating 64-gate case had an end-to-end ratio of about 1.06 but a
steady-state ratio below one. These are marginal and not a robust region-wide
advantage.

At 50,000 terms, 15% occupancy, 8 alternating segments, and 16 gates/segment:

- VectorPauliSum: 0.701 s construction-inclusive median.
- DenseLCB: 0.759 s construction-inclusive median.
- DenseLCB speedup: 0.923, i.e. about 8% slower.
- DenseLCB allocation measurement: approximately 0.93 GB.
- Crossing fallbacks: zero.

The remaining 50k high-occupancy/high-gate tail and all 100k–1M cases are
listed in `resource_limits.jsonl` as resource-limited skips. No terms were
truncated to complete them.

## 4. Regrouping result

The corrected alternating experiment reports exactly `segments - 1` active
region changes, one regrouping per change, and zero crossing gates. This fixes
the central methodological error in the earlier experiment.

Regrouping time is separately reported from local-gate time. It increases
with operator size and is not eliminated by denser local blocks. No stable
alternating-support break-even point was observed in the measured matrix.
The favorable zero-regrouping cases are a best-case control, not evidence that
active-region switching is profitable.

## 5. Realistic intermediate operators

`realistic_intermediate.jsonl` captures exact intermediate operators produced
by upstream Ising and random circuits, converts those same operators into
DenseLCB, and applies the same continuation circuit.

Examples:

| Family/depth | Captured terms | Mean occupancy | VectorPS time | DenseLCB time | Conversion to LCB |
|---|---:|---:|---:|---:|---:|
| Ising/8 | 13,999 | 78.1% | 1.46 ms | 2.36 ms | 13.4 ms |
| Random/8 | 1,129 | 62.4% | 0.159 ms | 0.242 ms | 1.31 ms |
| Random/8 | 909 | 42.1% | 0.146 ms | 0.200 ms | 1.12 ms |

The high-occupancy Ising operator does not make DenseLCB faster on the
realistic continuation. Conversion cost is several times larger than the
continuation itself in these cases.

## 6. Complete-circuit control

The corrected complete-circuit control is in `complete_circuits.jsonl` and
covers Ising/random families, widths 5/8, depths 4/8, and two seeds. With a
fixed `[1,2]` region, 36 of 56 primitive gates per 5-qubit Ising depth-4
circuit crossed the region; for 8-qubit depth-8 Ising, 144 of 184 crossed it.
This is the expected unfavorable case for LCB. DenseLCB was slower in every
complete-circuit case: representative speedups (VectorPauliSum/DenseLCB)
were 0.12–0.13 for 5-qubit Ising depth 4, 0.04–0.08 for 8-qubit Ising depth
4–8, and approximately 0.21–0.35 for the random controls. The complete
control therefore confirms that a local-kernel result cannot be promoted to a
whole-circuit result when gate supports move across the fixed region.

## 7. Timing interpretation

V2 distinguishes:

1. Pure local-kernel runtime on preconstructed representations.
2. Regrouping-inclusive runtime on preconstructed representations.
3. Construction-inclusive end-to-end runtime.

The first category can show a DenseLCB microbenchmark win. The second shows
whether switching active regions is affordable. The third is the practical
workflow comparison. Claims in this report use the third category unless
explicitly labeled steady-state.

## 8. Limitations

- The main LCB region is `k=2`; larger dense blocks were not justified after
  the corrected k=2 results.
- SparseLCB uses sorted small vectors and remains slower because insertion and
  repeated local merging dominate.
- Dense block clearing and outside-group allocation can consume substantial
  memory; the 50k measurement approached the configured conservative budget.
- The 50k campaign is intentionally partial and the larger requested scales
  are explicit skips, not silently absent data.
- Construction cost depends on the practical source representation. V2 reports
  both conversion-inclusive and steady-state values so neither interpretation
  is hidden.

## 9. Final verdict

**DenseLCB provides only synthetic/microbenchmark advantages.** It can be
useful as an optional specialized kernel for long, stable-support runs at
moderate-to-high occupancy, but the advantage is fragile and disappears or
reverses when support changes, realistic intermediate operators, conversion,
or larger operators are included. No adaptive active-region policy should be
added on the basis of these results.

Previous results remain under `results/`; all corrected v2 data and plots are
under this directory.
