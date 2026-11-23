"""Generate standalone scientific SVG/PNG plots and report from measured CSV."""
from pathlib import Path
import csv, json, os, statistics, sys, platform
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / "work" / "matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAW = ROOT / "results" / "raw"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader((RAW / "summary.csv").open(newline="", encoding="utf-8")))
raw = list(csv.DictReader((RAW / "benchmarks.csv").open(newline="", encoding="utf-8")))
meta = json.loads((RAW / "metadata.json").read_text())
cfg = meta["config"]
(RAW / "plot_environment.json").write_text(json.dumps({"python":sys.version,"matplotlib":matplotlib.__version__,"platform":platform.platform()},indent=2)+"\n")
primary = [r for r in rows if r["suite"] == "primary"]
backends = ["PauliSum", "VectorPauliSum", "bestfirst"]
colors = dict(zip(backends, ["#1672a4", "#e09600", "#b33461"]))
key = lambda r: tuple(r[k] for k in ("nqubits", "family", "layers", "seed", "observable"))
num = lambda r, k: float(r[k]) if r.get(k) else float("nan")
keys = sorted({key(r) for r in primary}, key=lambda x:(int(x[0]), x[1], int(x[2]), int(x[3]), x[4]))
index = {k:i for i,k in enumerate(keys)}
plt.rcParams.update({"font.size": 9, "figure.dpi": 140, "axes.spines.top": False, "axes.spines.right": False})

def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / (name + ".svg"))
    fig.savefig(FIG / (name + ".png"))
    plt.close(fig)

def decorate(ax):
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([f"{k[0]} {k[1][0]} L{k[2]} s{k[3]} {k[4].split('_')[0]}" for k in keys], rotation=90, fontsize=6)
    ax.set_xlabel("qubits / family (i=Ising, r=random) / layers / seed / observable")
    ax.grid(axis="y", alpha=.2)

for metric, label, name in [("median_ms","Median simulation time (ms, log scale)","runtime"),
                            ("allocated_bytes","Median allocated bytes (log scale)","allocations"),
                            ("peak_rss_mib","Process lifetime peak RSS (MiB; includes JIT/warmup)","memory")]:
    fig, ax = plt.subplots(figsize=(16,5))
    for backend in backends:
        subset = [r for r in primary if r["backend"] == backend]
        exact = [r for r in subset if r["status"] == "exact"]
        incomplete = [r for r in subset if r["status"] != "exact"]
        ax.scatter([index[key(r)] for r in exact], [num(r,metric) for r in exact], s=18, color=colors[backend], label=backend)
        if incomplete:
            ax.scatter([index[key(r)] for r in incomplete], [num(r,metric) for r in incomplete], marker="x", s=40, color=colors[backend], label=backend+" incomplete")
    if metric != "peak_rss_mib": ax.set_yscale("log")
    ax.set_ylabel(label); decorate(ax); ax.legend()
    save(fig,name)

fig,ax=plt.subplots(figsize=(16,5))
for backend in backends:
    sub=[r for r in primary if r["backend"]==backend]
    metric="peak_frontier" if backend=="bestfirst" else "peak_terms"
    ax.scatter([index[key(r)] for r in sub],[max(1,num(r,metric)) for r in sub],color=colors[backend],s=20,label=backend+(" frontier paths" if backend=="bestfirst" else " merged terms"))
ax.set_yscale("log"); ax.set_ylabel("Peak count (different quantities; see report)"); decorate(ax); ax.legend()
save(fig,"growth")

fig,axes=plt.subplots(1,2,figsize=(12,5))
budget=[r for r in rows if r["suite"]=="budget"]
for n in sorted({r["nqubits"] for r in budget}, key=int):
    for fam in ("ising","random"):
        sub=sorted([r for r in budget if r["nqubits"]==n and r["family"]==fam],key=lambda r:num(r,"runtime_cap_seconds"))
        axes[0].plot([num(r,"runtime_cap_seconds")*1000 for r in sub],[max(1e-16,num(r,"absolute_error")) for r in sub],"o-",label=f"n={n} {fam}")
        axes[1].plot([num(r,"median_ms") for r in sub],[max(1e-16,num(r,"absolute_error")) for r in sub],"o-",label=f"n={n} {fam}")
for ax in axes:
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=.2); ax.set_ylabel("Median absolute error (plot floor 1e-16)")
axes[0].set_xlabel("Requested BFS runtime cap (ms)")
axes[1].set_xlabel("Measured median runtime (ms)"); axes[1].legend(fontsize=8)
save(fig,"accuracy_budget")

ratios=[]
fig,ax=plt.subplots(figsize=(16,5))
for backend in backends[:2]:
    xs=[]; ys=[]
    for k in keys:
        bf=next((r for r in primary if key(r)==k and r["backend"]=="bestfirst" and r["status"]=="exact"),None)
        up=next((r for r in primary if key(r)==k and r["backend"]==backend and r["status"]=="exact"),None)
        if bf and up:
            ratio=num(up,"median_ms")/num(bf,"median_ms")
            xs.append(index[k]); ys.append(ratio); ratios.append((backend,k,ratio))
    ax.scatter(xs,ys,s=20,color=colors[backend],label=backend+" / bestfirst")
ax.axhline(1,color="black",linestyle="--",linewidth=1); ax.set_yscale("log")
ax.set_ylabel("Upstream time / BFS time; above 1 favors BFS; exact tasks only")
decorate(ax); ax.legend(); save(fig,"performance_ratio")

def fmt(x):
    if not x: return "—"
    try: return f"{float(x):.4g}"
    except ValueError: return str(x)

lines=["# Best-first Pauli propagation: measured report", "",
       f"Generated from {len(raw)} repeat/failure records. Primary timings use {cfg['repeats']} repetitions after a full configuration warmup. No compilation, circuit conversion, reference calculation or plotting is timed.", "",
       f"Environment: Julia {meta['julia_version']}, {meta['os']} ({meta['machine']}), {meta['cpu']}; {meta['logical_cpu_count']} logical CPUs detected, {meta['total_memory_bytes']/2**30:.2f} GiB RAM. One Julia thread and one BLAS thread. UTC run timestamp: {meta['timestamp']}.", "",
       f"OS detail: {meta.get('os_version','unavailable')}. Detected GPU (unused): {meta.get('gpu','unavailable')}.", "",
       "Upstream is registered PauliPropagation 0.9.2, pinned by Manifest.toml. Full metadata, source verification, parameter vectors, explicit circuits, repeat-level data and all package versions are under results/raw and Manifest.toml.", "",
       "## Observations", ""]
for backend in backends[:2]:
    rr=[r[2] for r in ratios if r[0]==backend]
    if rr:
        lines.append(f"- Against {backend}: BFS is faster in {sum(r>1 for r in rr)}/{len(rr)} completed equivalent cases. Median upstream/BFS runtime ratio is {statistics.median(rr):.3g}; range {min(rr):.3g}–{max(rr):.3g}. Ratios below 1 mean BFS is slower.")
incomplete=[r for r in primary if r["status"]!="exact"]
incomplete_reasons=Counter(r['reason'] for r in raw if r['suite']=='primary' and r['status']!='exact')
lines += [f"- {len(incomplete)} primary configurations include incomplete results. These are excluded from performance ratios and remain in every results table.",
          f"- Incomplete primary repeat reasons: {dict(incomplete_reasons)}. All upstream primary runs completed exactly in this recorded suite." if all(r['status']=='exact' for r in primary if r['backend']!='bestfirst') else f"- Incomplete primary repeat reasons: {dict(incomplete_reasons)}.",
          f"- Largest absolute error among validated exact repeats: {max((num(r,'absolute_error') for r in raw if r['status']=='exact' and r['absolute_error']),default=float('nan')):.3g}.",
          f"- Peak observed BFS frontier: {max((num(r,'peak_frontier') for r in primary if r['backend']=='bestfirst'),default=0):.0f} paths. Each live path occupies 48 payload bytes; vector capacity, allocator overhead and runtime memory are additional.", "",
          "## Interpretation and limits", "",
          "BFS explores unmerged paths. Upstream combines equal Pauli terms after each gate, allowing cancellation and reconvergence. BFS pays heap insertion/removal, bound bookkeeping and per-operation time checks, and revisits commuting gates separately along every path. Tiny sparse cases can favor the small BFS kernel; larger branching cases can strongly favor merging. These measurements describe this pilot grid, not general 20-qubit scalability.", "",
          "Several shallow timings are only a few microseconds; relative ratios there are sensitive to scheduling and timer overhead. Worker/backend order was fixed, not randomized. Width alone does not rank difficulty: a local observable has a depth-limited light cone, random draws vary with width, and gate/observable placement changes reconvergence. A 20-qubit case can therefore be easier than a 10-qubit case at the same layer count.", "",
          "The implemented priority is absolute path coefficient. It bounds the net remaining contribution of that subtree at fixed angles but does not predict its value. The sum of unfinished absolute coefficients gives an absolute error bound (with floating-point allowance); this is not a parameter-uniform LOWESA guarantee. Accuracy need not improve monotonically as budgets grow because partial contributions can cancel.", "",
          f"Primary results have no coefficient or weight truncation. BFS caps are {cfg['max_expansions']:,} expansions, {cfg['max_frontier_nodes']:,} frontier paths and {cfg['max_runtime_seconds']:g} seconds. These limits are not equivalent to upstream term counts. Every worker additionally has a {cfg['worker_rss_limit_mib']/1024:g} GiB sampled RSS cap and a {cfg['worker_case_timeout_seconds']:g}-second case timeout including warmup/diagnostics. The supervisor runs one worker at a time and polls every 50 ms; caps can overshoot between samples. External failures are recorded and subsequent cases continue.", "",
          f"The budget figure is a BFS-only anytime experiment at {', '.join(f'{b*1000:g}' for b in cfg['budget_seconds'])} ms with the same node caps. It is not an equal-cost comparison against upstream; upstream has no matching partial-tree result. No speedup is claimed between partial and exact runs.", "",
          "RSS is the process-lifetime high-water from Sys.maxrss, including Julia, compilation, warmup and earlier configurations in the same backend/size worker. It is not per-call incremental memory and must not be interpreted as algorithm-only memory. Allocated bytes measure allocation traffic, not peak live memory. External RSS sampling provides an independent safety guard. Peak upstream terms are collected in a separate untimed @peakpaulis pass after gate merging; BFS frontier counts represent unfinished paths. Upstream path/leaf/queue metrics are missing because they have no corresponding meaning.", "",
          "Exact 5-qubit references use independent dense state-vector evolution and agreement with both upstream backends. Larger references require agreement of both exact upstream backends; these share library kernels and are not fully independent validations. Every exact output and interrupted bound is checked against its available reference at atol=1e-10. Values at the error plot floor are displayed at 1e-16; CSV retains actual errors.", "",
          "Mean frontier occupancy, heap push/pop counts, allocations and GC time are recorded. Median occupancy and isolated heap-operation time are not measured; instrumenting each heap operation would change the timings. Final unique terms are available for upstream; BFS deliberately does not construct a merged final operator. This is an expectation-value benchmark, so that stronger upstream output is a relevant asymmetry.", "",
          "Future work: profile heap costs separately; evaluate priorities that account for remaining gates; investigate cancellation-aware merging as a separately specified algorithm; use repeated independent sessions to assess thermal/scheduling variation. No GPU, learned priorities or hybrid DFS were added.", "",
          "## Figures", ""]
profile_path=RAW/'profile_search.json'
if profile_path.exists():
    prof=json.loads(profile_path.read_text()); counts=prof['counts']
    note=f"A separate warmed sampling profile of {prof['circuit']} found heap frames in {counts['heap_samples']}/{counts['samples']} samples ({counts['heap_samples']/max(1,counts['samples']):.1%}). This supports queue overhead as a bottleneck on that workload; inclusive stack samples overlap and are not a precise fraction of wall-clock time. See raw/profile_search.txt and raw/profile_search.json."
    position=lines.index('## Figures')
    lines[position:position]=[note,'']
for name in ("runtime","allocations","memory","growth","accuracy_budget","performance_ratio"):
    lines += [f"![{name}](figures/{name}.svg)", ""]
for n in sorted({r['nqubits'] for r in primary},key=int):
    lines += [f"## {n} qubits", "", "| Family | Layers | Seed | Observable | Backend | Status | Median ms | IQR ms | RSS MiB† | Alloc bytes | Abs error | Bound | Peak paths/terms |", "|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in primary:
        if r['nqubits']!=n: continue
        fields=[r[k] for k in ('family','layers','seed','observable','backend','status')]+[fmt(r.get(k,'')) for k in ('median_ms','iqr_ms','peak_rss_mib','allocated_bytes','absolute_error','frontier_bound')]+[fmt(r['peak_frontier'] or r['peak_terms'])]
        lines.append('| '+' | '.join(fields)+' |')
    lines += ["", "† Process-lifetime high-water, not per-call incremental RSS.", ""]
(ROOT / 'results' / 'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
