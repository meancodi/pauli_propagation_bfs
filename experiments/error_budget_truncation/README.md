# Error-budgeted coefficient truncation

This is an isolated experiment for the research question: does online allocation of a global discarded-coefficient-mass budget improve runtime or memory at comparable certified error?

The certified implementation owns the full propagation loop. It merges equal Pauli words before every truncation checkpoint and records the exact discarded mass at that checkpoint. Upstream `min_abs_coeff` is not used by the certified policies, so there is no untracked implicit truncation. Policy D (`cost_aware`) uses only the current retained term count and recent term-count growth; it has no lookahead or preprocessing pass.

## Reproduce

From the repository root:

```powershell
julia --project=experiments/error_budget_truncation -e "using Pkg; Pkg.instantiate()"
julia --project=experiments/error_budget_truncation experiments/error_budget_truncation/test/runtests.jl
julia --project=experiments/error_budget_truncation --threads=1 experiments/error_budget_truncation/scripts/benchmark.jl --smoke
julia --project=experiments/error_budget_truncation experiments/error_budget_truncation/scripts/summarize.jl
python experiments/error_budget_truncation/scripts/plot_results.py
```

Append the bounded representative 20-qubit run after the smoke run:

```powershell
julia --project=experiments/error_budget_truncation --threads=1 experiments/error_budget_truncation/scripts/benchmark.jl --representative20 --append
```

Omit `--smoke` for the configured 5/10/20-qubit sweep. The benchmark writes append-only JSONL records to `results/raw/benchmark.jsonl`; incomplete records are retained and never treated as complete certificates. `epsilon_total=0` is the exact control. Fixed thresholds are a separate trade-off curve and have no certified global budget unless independently compared by observed error.

## Accounting definition

At checkpoint k, after merging, the discarded operator is the sum of discarded Pauli terms and the ledger adds `delta_k = sum(abs(coefficient))`. The reported certificate is `sum(delta_k)`, with conservative floating-point slack only for comparison. Unused allowance carries forward. A resource-limited run is marked `incomplete`.

## Scope and decision gate

Policies A--C are the required simple baselines. Policy D is intentionally transparent and is only meaningful if it wins repeatedly at the same certified budget without selection overhead erasing the gain. The report should preserve negative results.
