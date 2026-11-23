# Hierarchical merge-and-truncate v2 plan

The v2 study preserves the existing smoke results and adds a separate result
namespace under `results/hierarchical_v2/`.

1. Keep the existing `HierarchicalPauli` API stable and add a v2 module with
   explicit execution-plan nodes. A composed block is either materialized or a
   sequential plan; resource rejection never rebuilds a full exact map.
2. Add exact local maps with uniform telemetry: generated contributions,
   duplicate aggregation, true near-zero cancellation, peak temporary terms,
   peak map entries, and cumulative tree costs.
3. Add deterministic column truncation at internal nodes with modes `none`,
   `root_only`, `all_internal`, `alternate_levels`, `size_triggered`, and
   `selected_levels`, plus coefficient and per-column term limits.
4. Carry per-column infinity-norm error bounds through child composition. These
   are conservative bounds for the sparse-column approximation and are reported
   separately from observed expectation error.
5. Implement genuinely distinct left-deep, balanced, and support-constrained
   planning. Support-constrained planning uses only support and estimated
   column/storage costs; it never trial-composes candidates.
6. Add an on-demand column cache with explicit missing/unavailable states.
7. Add bounded resource estimates and a resumable, single-process progressive
   campaign. Every v2 record is append-only, keyed by a configuration id, and
   records completion, skip, timeout, or resource-limit status.
8. Run tests first, then moderate unattended tiers, then targeted stress cases;
   regenerate v2 figures and reports without overwriting earlier data.

The v2 result is scientific feasibility evidence. Any unavailable rigorous
bound or soft platform memory monitor is explicitly labelled rather than
presented as a guarantee.
