# LCB benchmark correction plan

The original results remain in `results/` and are not overwritten. Corrected
raw data will be written to `results/validated_v2/`.

## Defects being corrected

1. Synthetic masks were sampled from all 64 bits even when `n < 64`.
2. Group generation used an unbounded random loop and did not report the full
   occupancy distribution or feasibility bound.
3. The previous regrouping loop changed the active region but kept gates on
   qubits 1 and 2, so it measured crossing fallbacks rather than local-gate
   amortization.
4. Synthetic allocation timing ran on the operator after the timed sequence,
   not on an independent copy of the initial state.
5. The prior upstream circuit comparison did not provide a matched,
   preconstructed `VectorPauliSum` for the alternating-support experiment.
6. The previous runtime plot connected unrelated width/depth/family points.

## Corrected protocol

The v2 generator uses a bounded base-four enumeration of outside configurations
under the declared n-qubit mask. Break-even cases use alternating `[1,2]` and
`[3,4]` segments whose gate supports match their active region; a zero-switch
control keeps all gates in `[1,2]`. Each matched case has independent copies,
warm-up, repeated measurements, construction-inclusive timing, and paired
final-operator/expectation validation. Upstream propagation uses the public
`PP.propagate(...; min_abs_coeff=0.0, max_weight=Inf, thread=false)` path with
an actual `VectorPauliSum` object and Float64 coefficients.

The campaign is progressive and records skipped resource-limited configurations
explicitly. No truncation or adaptive policy is introduced.
