# Favorable structured-circuit study

This isolated study preserves the earlier hierarchical experiments and tests
five deliberately structured families: near-inverse echoes, repeated local
blocks, Clifford sandwiches, perturbed Trotter blocks, and periodic finite-order
motifs.

Run the reproducible checks and screen:

```powershell
julia --project=experiments/hierarchical_pauli experiments/hierarchical_favorable_circuits/test/runtests.jl
julia --project=experiments/hierarchical_pauli --threads=1 experiments/hierarchical_favorable_circuits/benchmarks/screen.jl
python experiments/hierarchical_favorable_circuits/benchmarks/generate_reports.py
```

The screen writes `results/favorable_v1/screen_v2.jsonl`. The curated moderate
runner is resumable and uses deterministic SHA-1 configuration IDs:

```powershell
julia --project=experiments/hierarchical_pauli --threads=1 experiments/hierarchical_favorable_circuits/benchmarks/run_campaign.jl
```

It writes durable `campaign.jsonl` records and progress logs. A campaign lock
prevents overlapping writers; use `--force-lock` only after verifying that no
campaign process remains.

The screen and partial campaign found no reproducible recursive advantage over
both optimized propagation and simpler block/caching alternatives. See
[`FINAL_REPORT.md`](FINAL_REPORT.md).
