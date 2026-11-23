# Final report

The measured report is maintained at
[`results/favorable_v1/FINAL_REPORT.md`](results/favorable_v1/FINAL_REPORT.md).

The five-family structural screen passed dense-reference correctness checks but
did not identify a family where recursive merge-and-truncate beat both
block-boundary propagation and the optimized upstream backend after construction
cost was included. Repeated local blocks remain the most plausible restricted
follow-up because they can amortize exact cached maps; the current screen did
not yet demonstrate that advantage for recursive map powers.
