# Best-first Pauli propagation

Independent Julia best-first search for fixed-angle Pauli expectation values,
with tests and reproducible comparisons against official PauliPropagation.jl.

See the [experiment README](experiments/best_first_pauli/README.md) for setup,
API, benchmark commands and measurement definitions, and the
[measured report](experiments/best_first_pauli/results/report.md) for results.

The DFS equal-budget extension and its independent outputs are documented in the
[experiment README](experiments/best_first_pauli/README.md) and
[DFS comparison report](experiments/best_first_pauli/results_dfs_vs_bestfirst/report.md).

The original [design document](best_first_pauli_propagation_design.md) is preserved.

The error-budget truncation experiment is in [experiments/error_budget_truncation](experiments/error_budget_truncation/), with reproducible tests, benchmark scripts, raw JSONL data, figures, and a preliminary report.
