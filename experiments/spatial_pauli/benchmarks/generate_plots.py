"""Generate only charts supported by measured CSV columns."""
import csv
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..", "results")
paths = sorted(glob.glob(os.path.join(ROOT, "raw", "benchmark_*.csv")))
if not paths:
    raise SystemExit("No benchmark CSV found")
merged = {}
for path in paths:
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if not r.get("CPU"):
                continue
            key = (r["family"], r["n"], r["depth"], r["observable"], r["coupling"], r["method"])
            merged[key] = r
data = [r for r in merged.values() if r["status"] == "ok"]
figdir = os.path.join(ROOT, "figures")
os.makedirs(figdir, exist_ok=True)
colors = {"S0": "#333333", "S1": "#1f77b4", "S2": "#ff7f0e"}
families = sorted({r["family"] for r in data})
ns = sorted({int(r["n"]) for r in data})

def selected():
    return [r for r in data if r["observable"] == "local_sum" and r["depth"] == "1"
            and int(r["coupling"]) == (1 if r["family"] in ("weak", "cross") else 0)]

for field, ylabel, name in [
    ("median_total_s", "Median total time (s)", "runtime.png"),
    ("peak_live_bytes", "Peak live representation estimate (bytes)", "memory.png"),
]:
    rows = selected()
    fig, axes = plt.subplots(1, len(families), figsize=(4 * len(families), 3.5), sharey=True)
    for ax, family in zip(axes, families):
        for method, offset in (("S0", -0.24), ("S1", 0), ("S2", 0.24)):
            vals = {int(r["n"]): float(r[field]) for r in rows if r["family"] == family and r["method"] == method}
            ax.bar([i + offset for i in range(len(ns)) if ns[i] in vals], [vals[n] for n in ns if n in vals], width=0.23, color=colors[method], label=method)
        ax.set_xticks(range(len(ns)), ns)
        ax.set_title(family)
        ax.set_xlabel("Qubits")
        ax.set_yscale("log")
    axes[0].set_ylabel(ylabel)
    axes[-1].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, name), dpi=160)
    plt.close(fig)

base = {(r["family"], r["n"], r["depth"], r["observable"]): r for r in data if r["method"] == "S0"}
fig, ax = plt.subplots(figsize=(10, 4))
labels, vals, shades = [], [], []
for r in selected():
    if r["method"] == "S0":
        continue
    b = base.get((r["family"], r["n"], r["depth"], r["observable"]))
    if b:
        labels.append(f'{r["family"]}/{r["n"]}/{r["method"]}')
        vals.append(100 * (1 - float(r["median_total_s"]) / float(b["median_total_s"])))
        shades.append(colors[r["method"]])
ax.bar(range(len(vals)), vals, color=shades)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_yscale("symlog", linthresh=100)
ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=7)
ax.set_ylabel("Runtime improvement vs S0 (%)")
fig.tight_layout()
fig.savefig(os.path.join(figdir, "improvement.png"), dpi=160)
plt.close(fig)

rows = [r for r in data if r["family"] == "ising" and r["n"] == "8" and r["observable"] == "local_sum" and r["method"] != "S0"]
if rows:
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for method in ("S1", "S2"):
        rr = sorted((r for r in rows if r["method"] == method), key=lambda r: int(r["depth"]))
        ax.plot([int(r["depth"]) for r in rr], [float(r["peak_components"]) for r in rr], marker="o", label=method)
    ax.set(xlabel="Circuit depth", ylabel="Maximum factorized components", title="8-qubit Ising, local-sum observable")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "rank_growth.png"), dpi=160)
    plt.close(fig)

# Coupling count is plotted only for a dedicated sweep at fixed n and depth.
coupling_rows = [r for r in data if r["family"] == "weak" and r["n"] == "8" and r["depth"] == "1" and r["observable"] == "local_sum"]
if len({r["coupling"] for r in coupling_rows}) >= 2:
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    for method in ("S0", "S1", "S2"):
        rr = sorted((r for r in coupling_rows if r["method"] == method), key=lambda r: int(r["coupling"]))
        axes[0].plot([int(r["coupling"]) for r in rr], [float(r["median_total_s"]) for r in rr], marker="o", label=method)
        axes[1].plot([int(r["coupling"]) for r in rr], [float(r["peak_live_bytes"]) for r in rr], marker="o", label=method)
    axes[0].set(xlabel="Cross-block interactions", ylabel="Median total time (s)")
    axes[1].set(xlabel="Cross-block interactions", ylabel="Peak live estimate (bytes)")
    axes[0].legend(); fig.tight_layout()
    fig.savefig(os.path.join(figdir, "coupling_sensitivity.png"), dpi=160)
    plt.close(fig)
print("Generated figures from", len(paths), "CSV files")
