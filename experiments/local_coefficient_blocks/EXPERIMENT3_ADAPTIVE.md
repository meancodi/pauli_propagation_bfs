# Experiment 3: adaptive active regions

Not implemented as a performance feature. The prerequisite fixed-region
end-to-end campaign did not show a competitive advantage over
`VectorPauliSum`, and regrouping is a full-term scan. Adding look-ahead would
therefore risk optimizing a synthetic local kernel while making complete
circuits slower.

The implementation deliberately retains the exact, configurable `regroup!`
primitive so a future experiment can add a policy after a stronger workload
with long stable gate-support runs is identified.
