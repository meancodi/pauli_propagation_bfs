# Favorable circuit families: measured report

The primary structural screen contains five families and paired backends with
three repetitions per observable. The moderate campaign was started with five
repetitions and durable JSONL output, but was intentionally stopped during the
expensive repeated-block tail after substantial echo/repeated records were
written. It is reported as partial, not as a completed Cartesian campaign.

## Screen numerical results

- `clifford_sandwich` / `B0_upstream`: 3 measurements; median runtime `0.0001674` s; median error `0.0`; nonzero-error cases `0`.
- `clifford_sandwich` / `B2_gatewise`: 3 measurements; median runtime `4.92e-05` s; median error `0.0`; nonzero-error cases `0`.
- `clifford_sandwich` / `B3_boundary`: 3 measurements; median runtime `5.21e-05` s; median error `0.0`; nonzero-error cases `0`.
- `clifford_sandwich` / `B4_cached_local`: 3 measurements; median runtime `0.002646` s; median error `0.0`; nonzero-error cases `0`.
- `clifford_sandwich` / `B5_balanced_hier`: 3 measurements; median runtime `0.0188644` s; median error `0.0`; nonzero-error cases `0`.
- `clifford_sandwich` / `B6_support_hier`: 3 measurements; median runtime `0.0045589` s; median error `0.0`; nonzero-error cases `0`.
- `clifford_sandwich` / `B8_ondemand`: 3 measurements; median runtime `0.0003512` s; median error `0.0`; nonzero-error cases `0`.
- `echo` / `B0_upstream`: 3 measurements; median runtime `0.0001718` s; median error `0.0`; nonzero-error cases `0`.
- `echo` / `B2_gatewise`: 3 measurements; median runtime `6.19e-05` s; median error `6.683373247873289e-18`; nonzero-error cases `0`.
- `echo` / `B3_boundary`: 3 measurements; median runtime `7.7e-05` s; median error `6.683373247873289e-18`; nonzero-error cases `0`.
- `echo` / `B4_cached_local`: 3 measurements; median runtime `0.0012873` s; median error `0.0`; nonzero-error cases `0`.
- `echo` / `B5_balanced_hier`: 3 measurements; median runtime `0.121824` s; median error `6.683373247873289e-18`; nonzero-error cases `0`.
- `echo` / `B6_support_hier`: 3 measurements; median runtime `0.0065235` s; median error `0.0`; nonzero-error cases `0`.
- `echo` / `B8_ondemand`: 3 measurements; median runtime `0.0006585` s; median error `0.0`; nonzero-error cases `0`.
- `periodic` / `B0_upstream`: 3 measurements; median runtime `0.0001333` s; median error `0.0`; nonzero-error cases `0`.
- `periodic` / `B2_gatewise`: 3 measurements; median runtime `3.29e-05` s; median error `0.0`; nonzero-error cases `0`.
- `periodic` / `B3_boundary`: 3 measurements; median runtime `1.91e-05` s; median error `0.0`; nonzero-error cases `0`.
- `periodic` / `B4_cached_local`: 3 measurements; median runtime `9.37e-05` s; median error `0.0`; nonzero-error cases `0`.
- `periodic` / `B5_balanced_hier`: 3 measurements; median runtime `0.0016047` s; median error `0.0`; nonzero-error cases `0`.
- `periodic` / `B6_support_hier`: 3 measurements; median runtime `0.0012326` s; median error `0.0`; nonzero-error cases `0`.
- `periodic` / `B8_ondemand`: 3 measurements; median runtime `0.0001975` s; median error `0.0`; nonzero-error cases `0`.
- `repeated` / `B0_upstream`: 3 measurements; median runtime `0.0006299` s; median error `0.0`; nonzero-error cases `0`.
- `repeated` / `B2_gatewise`: 3 measurements; median runtime `0.0039697` s; median error `0.0011358311739555482`; nonzero-error cases `3`.
- `repeated` / `B3_boundary`: 3 measurements; median runtime `0.0038125` s; median error `0.00015053832936848255`; nonzero-error cases `3`.
- `repeated` / `B4_cached_local`: 3 measurements; median runtime `0.0050044` s; median error `1.3877787807814457e-16`; nonzero-error cases `0`.
- `repeated` / `B5_balanced_hier`: 3 measurements; median runtime `0.557946` s; median error `0.0012126610583132191`; nonzero-error cases `3`.
- `repeated` / `B6_support_hier`: 3 measurements; median runtime `0.589157` s; median error `0.0012126610583132191`; nonzero-error cases `3`.
- `repeated` / `B7_repeated_powers`: 3 measurements; median runtime `0.278811` s; median error `3.469446951953614e-16`; nonzero-error cases `0`.
- `repeated` / `B8_ondemand`: 3 measurements; median runtime `0.0050696` s; median error `0.0`; nonzero-error cases `0`.
- `trotter` / `B0_upstream`: 3 measurements; median runtime `0.0005964` s; median error `0.0`; nonzero-error cases `0`.
- `trotter` / `B2_gatewise`: 3 measurements; median runtime `0.0001892` s; median error `0.0020863473156659193`; nonzero-error cases `3`.
- `trotter` / `B3_boundary`: 3 measurements; median runtime `0.0001675` s; median error `0.0020864340159241768`; nonzero-error cases `3`.
- `trotter` / `B4_cached_local`: 3 measurements; median runtime `0.0457303` s; median error `1.3877787807814457e-17`; nonzero-error cases `0`.
- `trotter` / `B5_balanced_hier`: 3 measurements; median runtime `0.0434347` s; median error `0.0010742853625476737`; nonzero-error cases `3`.
- `trotter` / `B6_support_hier`: 3 measurements; median runtime `0.0225671` s; median error `0.0010742503550259727`; nonzero-error cases `3`.
- `trotter` / `B8_ondemand`: 3 measurements; median runtime `0.0209589` s; median error `6.938893903907228e-18`; nonzero-error cases `0`.

## Formulation and correctness

Each local map stores sparse columns of `H_B(P)=B†PB`. A parent composes the
right-child column through the left child, merges identical Pauli keys, and may
truncate before its ancestor is constructed. B5/B6 use genuine internal
truncation; B4 caches exact local maps; B7 constructs repeated powers where
applicable; B8 uses explicit on-demand fallback leaves.

All five generators pass dense state-vector versus upstream `VectorPauliSum`
checks at small widths. Echo inverse and periodic controls pass their algebraic
tests, and the previous v1/v2 suites remain green.

## Final evaluation

The screen shows no family where B5 or B6 beats both B3 and optimized B0 after
map preparation is included. Repeated blocks are the most plausible caching
case, but B7/power preparation still carries overhead. Echo, Clifford-sandwich,
and exact periodic cases have simpler algebraic alternatives. Trotter circuits
show nonzero truncation sensitivity but not a recursive speed advantage.

## Recommendation

Do not adopt general recursive transfer-map construction by default. Restrict
further work to repeated local blocks with very high reuse, and compare against
exact cached single-block application and compiler simplification first.
Ordinary block-boundary truncation remains the more credible practical path.
The expensive moderate tail was stopped and is explicitly reported as partial;
no missing results are fabricated.
