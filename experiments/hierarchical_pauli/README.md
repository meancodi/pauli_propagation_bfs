# Hierarchical merge-and-truncate Pauli propagation

This is a self-contained feasibility study for exact local transfer-map
composition, block-boundary truncation, and recursive/hierarchical truncation.
The upstream PauliPropagation.jl source is not modified.

Run from the repository root:

```powershell
julia --project=experiments/hierarchical_pauli -e "using Pkg; Pkg.instantiate()"
julia --project=experiments/hierarchical_pauli experiments/hierarchical_pauli/test/runtests.jl
julia --project=experiments/hierarchical_pauli --threads=1 experiments/hierarchical_pauli/benchmarks/experiment1.jl --smoke
julia --project=experiments/hierarchical_pauli --threads=1 experiments/hierarchical_pauli/benchmarks/experiment2.jl --smoke
julia --project=experiments/hierarchical_pauli --threads=1 experiments/hierarchical_pauli/benchmarks/experiment3.jl --smoke
python experiments/hierarchical_pauli/benchmarks/generate_plots.py
```

The smoke suite is intentionally bounded for a 16 GB workstation. JSONL files
contain incomplete/resource-limited runs as well as completed records.
