# Local Coefficient Blocks: final report

## Verdict

LCB is useful only as an optional specialized kernel for repeated gates inside
a stable two-qubit region and moderate local occupancy. It is not a
general-purpose replacement for `PauliSum` or optimized `VectorPauliSum` on the
tested complete circuits.

## Answers to the research questions

1. Grouping can reduce local update work, but only when gates stay inside the
   region. Crossing/random gates pay a full regrouping/materialization cost.
2. In this `k=2` implementation, dense blocks were favorable around the
   measured 15% occupancy point; sparse blocks never beat dense blocks in the
   measured campaign. The crossover is workload and implementation dependent,
   not a universal occupancy constant.
3. Regrouping needs a long stable run to amortize its full scan. The completed
   campaign shows substantial cost even with 1–8 local-gate groups per segment;
   no defensible universal break-even count can be claimed from the resource-
   capped tail.
4. LCB is not competitive when support changes frequently.
5. An adaptive policy was not justified: switching cost dominates the observed
   end-to-end opportunity.
6. It did not outperform `VectorPauliSum` in the 5/8-qubit full-circuit
   campaign, despite a favorable local microbenchmark.
7. The design differs from the repository's hierarchical transfer-map
   experiment and upstream `TransferMap`/`TransferMapGate`: LCB stores an
   evolving sum partitioned by outside configuration, whereas transfer maps
   store reusable local input-to-output columns. Both exploit locality and
   neither should be treated as automatically novel. Existing vector-packed
   and transfer-map paths remain stronger baselines.

## Limitations

The current experiment uses `k=2`, Float64, one thread, and an exact
correctness-first boundary path. Dense blocks allocate `4^k` values per group;
sparse blocks use sorted small vectors and can pay repeated insertion costs.
The campaign is progressive rather than a full million-term Cartesian sweep;
all durable raw records and the stopped resource-limited portion are retained.

## Reproduction

From the repository root:

```text
julia --project=experiments/hierarchical_pauli experiments/local_coefficient_blocks/test/runtests.jl
julia --project=experiments/hierarchical_pauli experiments/local_coefficient_blocks/benchmarks/synthetic.jl
julia --project=experiments/hierarchical_pauli experiments/local_coefficient_blocks/benchmarks/upstream.jl
julia --project=experiments/hierarchical_pauli experiments/local_coefficient_blocks/benchmarks/regrouping.jl
python experiments/local_coefficient_blocks/benchmarks/generate_plots.py
```

Existing repository tests and experiments were not modified.
