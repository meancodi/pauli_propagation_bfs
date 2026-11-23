# Experiment 2: active-region regrouping

`regroup!(x, S)` preserves every term and rebuilds the outside-key partition.
It is exact and supports noncontiguous regions. Boundary gate application uses
the same operation internally.

The campaign alternates `[1,2]` and `[3,4]` over eight segments and measures
1/2/4/8/16/32/64 local gates per region where resources permit. Completed
records are in `results/regrouping.jsonl`; the long 50k tail at 16–64 gates was
stopped after the completed 1/2/4/8 points to respect the resource budget.

At 50k terms, eight segments took 1.05 s dense versus 1.60 s sparse with one
local-gate group per segment, and 8.51 s dense versus 12.49 s sparse at eight
groups. Regrouping remains a full operator scan, so increasing the number of
local gates is the only available amortization mechanism. The data supports
using LCB only for long runs of gates with stable support; it does not support
switching regions after every gate.

No adaptive policy was enabled after this phase: the measured end-to-end
upstream advantage was absent, and a policy would add switching decisions on
top of an already expensive regrouping path.
