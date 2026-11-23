"""Audit coverage, repeats, reference errors and interruption bounds in saved data."""
from pathlib import Path
from collections import defaultdict
import csv
import json
import itertools
import math

root = Path(__file__).resolve().parents[1]
raw = root / "results" / "raw"
config = json.loads((raw / "metadata.json").read_text())["config"]
rows = list(csv.DictReader((raw / "benchmarks.csv").open(newline="", encoding="utf-8")))
assert len({r["run_id"] for r in rows}) == len(rows), "Duplicate run IDs"
groups = defaultdict(list)
for row in rows:
    if row["suite"] == "primary":
        k = (int(row["nqubits"]), row["family"], int(row["layers"]), int(row["seed"]), row["observable"], row["backend"])
        groups[k].append(row)
    if row["elapsed_ms"]:
        assert float(row["elapsed_ms"]) > 0
    if row["absolute_error"]:
        error = float(row["absolute_error"])
        assert math.isfinite(error)
        if row["status"] == "exact":
            assert error <= 1e-10, row["run_id"]
        elif row["backend"] == "bestfirst" and row["frontier_bound"]:
            assert error <= float(row["frontier_bound"]) + 1e-10, row["run_id"]
    if row["backend"] == "bestfirst" and row["peak_frontier"]:
        assert int(row["peak_frontier"]) <= int(row["max_frontier_nodes"])
        assert int(row["expanded_nodes"]) <= int(row["max_expansions"])
        if row["status"] == "exact":
            assert float(row["frontier_bound"]) == 0
            assert int(row["frontier_nodes"]) == 0
    if row["threads"]:
        assert row["threads"] == "1"
for k in itertools.product(config["nqubits"], config["families"], config["layers"],
                           config["seeds"], config["observables"], ("PauliSum", "VectorPauliSum", "bestfirst")):
    assert k in groups, f"Missing configuration: {k}"
    group = groups[k]
    external = any(r["reason"].startswith("external_") for r in group)
    if not external:
        assert len(group) == config["repeats"], f"Repeat count: {k}"
        assert all(r["reference"] for r in group), f"Missing reference: {k}"
assert all((root / "results" / "figures" / (name + ".svg")).is_file()
           for name in ("runtime", "memory", "allocations", "growth", "accuracy_budget", "performance_ratio"))
print(f"Validated {len(rows)} records across {len(groups)} primary backend configurations; coverage, repeats, errors and bounds passed.")
