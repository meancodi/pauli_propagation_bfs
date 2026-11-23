# Corrected v2 benchmark report

Campaign data: `D:\GitHub_Repo\pauli_propagation_bfs\experiments\hierarchical_pauli\results\hierarchical_v2\benchmark.jsonl`. The runner used one Julia thread and one BLAS thread, durable JSONL output, resume-by-config-id, and an explicit warning that Windows enforcement is admission/Julia-limit based rather than a hard RSS kill.

B0: 1408 completed records, median runtime 7.53e-05s, median observed error 0.0, nonzero-error cases 0, internal discarded-column events 0, resource events 0.
B1: 1408 completed records, median runtime 7.245e-05s, median observed error 0.0, nonzero-error cases 456, internal discarded-column events 0, resource events 0.
B2: 1408 completed records, median runtime 5.975e-05s, median observed error 0.0, nonzero-error cases 456, internal discarded-column events 0, resource events 0.
B3: 1408 completed records, median runtime 0.0070086s, median observed error 0.0, nonzero-error cases 0, internal discarded-column events 0, resource events 20000.
B4: 1408 completed records, median runtime 0.00758365s, median observed error 0.0, nonzero-error cases 190, internal discarded-column events 4340096, resource events 20000.
B5: 1408 completed records, median runtime 0.0022237s, median observed error 0.0, nonzero-error cases 6, internal discarded-column events 701184, resource events 79504.
B6: 1408 completed records, median runtime 0.0028737s, median observed error 0.0, nonzero-error cases 145, internal discarded-column events 1431936, resource events 31360.
B7: 1408 completed records, median runtime 0.0002323s, median observed error 0.0, nonzero-error cases 0, internal discarded-column events 0, resource events 173088.

## Experiment 2
B1 is gate-wise truncation and B2 is ordinary four-gate boundary truncation. Curves are grouped by family and width.
