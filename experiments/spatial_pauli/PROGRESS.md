# Progress

- Inspected the generic `PauliWord`, rotation, Clifford, circuit, and reference kernels in `best_first_pauli`.
- Implemented exact two-partition and tree leaf factorizations with parity-sector boundary rotations and letter-sector CNOTs.
- Added exact aggregation, resource limits, and coefficient-level tests.
- Approximate compression remains disabled pending evidence of a useful exact factorization advantage.
- Final validation: 157/157 tests passed, including upstream and dense expectation checks.
- Corrected benchmark: 333 distinct method records, 222 paired comparisons, five generated figures. Only 4 paired spatial runs beat S0 in runtime; a repeated 8-qubit cross-coupled ZZ case is the narrow positive result. See `results/REPORT.md`.
- Exploratory CSVs without a `CPU` column used timing contaminated by per-gate memory inspection. They are preserved but excluded from the corrected report and plots.
