"""Plot only exact, matched V2 records; never read V1 CSVs."""
import csv
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "results", "v2"))
raw = os.path.join(ROOT, "raw")
figdir = os.path.join(ROOT, "figures")
os.makedirs(figdir, exist_ok=True)
paths = []
for mode in ("pilot", "targeted", "seeded", "final", "coupling"):
    files = sorted(glob.glob(os.path.join(raw, f"v2_{mode}_*.csv")))
    if files:
        paths.append(files[-1])
rows = []
for path in paths:
    with open(path, newline="") as f:
        rows.extend(csv.DictReader(f))
if not rows:
    raise SystemExit("No V2 data")

CONFIG = ("version", "family", "n", "layers", "gate_count", "observable", "coupling",
          "periodic", "seed", "fingerprint", "partition_scheme", "leaf_size", "method")
dedup = {tuple(r[k] for k in CONFIG): r for r in rows}
rows = list(dedup.values())
good = [r for r in rows if r["status"] == "completed" and float(r["expectation_error"]) <= 1e-8]
PAIR = ("version", "family", "n", "layers", "gate_count", "observable", "coupling",
        "periodic", "seed", "fingerprint")
groups = defaultdict(dict)
for r in good:
    groups[tuple(r[k] for k in PAIR)][r["method"]] = r
paired = [(k, v) for k, v in groups.items() if "S0" in v]
methods = ("S0", "S0_optimized", "S1", "S2", "S3")
colors = dict(zip(methods, ("#333333", "#777777", "#1f77b4", "#ff7f0e", "#2ca02c")))

def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, name), dpi=170)
    plt.close(fig)

def grouped(field, ylabel, filename):
    selected = [r for r in good if r["n"] == "8" and r["observable"] == "factorized" and r["layers"] == "1"]
    families = sorted({r["family"] for r in selected})
    if not families:
        return
    fig, ax = plt.subplots(figsize=(max(8, len(families) * 1.65), 4))
    for mi, method in enumerate(methods):
        values = []
        for family in families:
            rr = [float(r[field]) for r in selected if r["family"] == family and r["method"] == method and r[field]]
            values.append(sorted(rr)[len(rr)//2] if rr else float("nan"))
        ax.bar([i + (mi-2)*0.16 for i in range(len(families))], values, width=0.15, label=method, color=colors[method])
    ax.set_xticks(range(len(families)), families, rotation=25)
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    ax.legend(ncol=5, fontsize=8)
    save(fig, filename)

grouped("median_total_s", "Median total runtime (s)", "runtime.png")
grouped("peak_retained_bytes", "Peak retained Julia object size (bytes)", "memory.png")
grouped("allocated_bytes", "Maximum allocated bytes across five timed runs", "allocations.png")

improvements = []
for key, v in paired:
    for method in methods[1:]:
        if method in v:
            improvements.append((f'{key[1]}/{key[2]}/{method}',
                                 100*(1-float(v[method]["median_total_s"])/float(v["S0"]["median_total_s"])),method))
if improvements:
    # One representative seed per family/method keeps labels readable.
    chosen = {}
    for label, value, method in improvements:
        family = label.split("/")[0]
        chosen.setdefault((family, method), (label, value, method))
    fig, ax = plt.subplots(figsize=(10, 4))
    items = list(chosen.values())
    ax.bar(range(len(items)), [x[1] for x in items], color=[colors[x[2]] for x in items])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_yscale("symlog", linthresh=100)
    ax.set_xticks(range(len(items)), [x[0] for x in items], rotation=90, fontsize=7)
    ax.set_ylabel("Runtime improvement vs matched S0 (%)")
    save(fig, "improvement.png")

fig, ax = plt.subplots(figsize=(6, 4))
for method in ("S1", "S2"):
    rr = [r for r in good if r["n"] == "8" and r["observable"] == "factorized" and r["family"] == "cross" and r["method"] == method]
    depths = sorted({int(r["layers"]) for r in rr})
    ys = [max(float(r["peak_components"]) for r in rr if int(r["layers"]) == d) for d in depths]
    if depths:
        ax.plot(depths, ys, marker="o", label=method)
ax.set(xlabel="Circuit layers", ylabel="Peak component count", title="Cross-partition family, 8 qubits")
ax.legend()
save(fig, "component_growth.png")

coupling = [r for r in good if r["family"] == "weak" and r["n"] == "8" and r["layers"] == "1" and r["observable"] == "factorized"]
if len({r["coupling"] for r in coupling}) >= 2:
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    for method in methods:
        rr = sorted((r for r in coupling if r["method"] == method), key=lambda r: int(r["coupling"]))
        if rr:
            axes[0].plot([int(r["coupling"]) for r in rr], [float(r["median_total_s"]) for r in rr], marker="o", label=method)
            axes[1].plot([int(r["coupling"]) for r in rr], [float(r["peak_retained_bytes"]) for r in rr], marker="o", label=method)
    axes[0].set(xlabel="Boundary interactions", ylabel="Median total runtime (s)")
    axes[1].set(xlabel="Boundary interactions", ylabel="Peak retained bytes")
    axes[0].legend(fontsize=7)
    save(fig, "coupling_sensitivity.png")

points = [(float(v["S2"]["median_total_s"]),float(v["S3"]["median_total_s"]))
          for _,v in paired if "S2" in v and "S3" in v]
if points:
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter([p[0] for p in points],[p[1] for p in points],s=12,alpha=0.7)
    lo=min(min(p) for p in points); hi=max(max(p) for p in points)
    ax.plot([lo,hi],[lo,hi],color="black",linewidth=0.8)
    ax.set(xlabel="S2 median total runtime (s)",ylabel="S3 median total runtime (s)")
    ax.set_xscale("log"); ax.set_yscale("log")
    save(fig,"s2_vs_s3.png")

fig, ax = plt.subplots(figsize=(6, 4))
for method in ("S0", "S1", "S2", "S3"):
    rr = [r for r in good if r["method"] == method and r["n"] == "8"]
    if rr:
        ax.scatter([float(r["median_init_s"]) for r in rr],
                   [float(r["median_gate_s"]) for r in rr],label=method,s=16,alpha=0.7)
ax.set(xlabel="Initialization time (s)",ylabel="Gate phase time (s)")
ax.set_xscale("symlog",linthresh=1e-6); ax.set_yscale("log")
ax.legend(fontsize=7)
save(fig,"construction_vs_propagation.png")
print(f"V2 figures: {len(good)} completed rows, {len(paired)} matched configurations")
