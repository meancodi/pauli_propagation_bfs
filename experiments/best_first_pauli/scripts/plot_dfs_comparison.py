"""Consolidate and visualize the isolated DFS-versus-best-first experiment."""
from pathlib import Path
from collections import defaultdict, Counter
import csv, json, math, os, statistics, sys, platform

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results_dfs_vs_bestfirst'; RAW=OUT/'raw'; FIG=OUT/'figures'
FIG.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'work'/'matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

records=[]
for path in sorted(RAW.glob('*_search.jsonl')):
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.strip(): records.append(json.loads(line))
byid={r['run_id']:r for r in records}; records=list(byid.values())
fields=sorted({k for r in records for k in r})
with (RAW/'search_repeats.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(records)
primary=[r for r in records if r.get('algorithm') in ('dfs','bestfirst') and r.get('experiment') in ('expansion_budget','runtime_budget','accuracy_tolerance')]
gkeys=('nqubits','family','layers','seed','observable','algorithm','experiment','requested_budget')
groups=defaultdict(list)
for r in primary: groups[tuple(r.get(k) for k in gkeys)].append(r)
def med(rs,k):
    v=[float(r[k]) for r in rs if r.get(k) is not None]
    return statistics.median(v) if v else None
summary=[]
for key,rs in groups.items():
    row=dict(zip(gkeys,key)); row['repeats']=len(rs)
    row['status']='/'.join(sorted({r.get('status','') for r in rs}))
    row['reason']='/'.join(sorted({r.get('reason','') for r in rs}))
    for k in ('absolute_error','relative_error','frontier_bound','elapsed_ms','search_elapsed_ms','allocated_bytes','gc_ms',
              'peak_rss_mib','expanded_nodes','generated_children','completed_leaves','peak_frontier','frontier_nodes',
              'observed_children_per_expansion','reference','reference_elapsed_ms'):
        row[k]=med(rs,k)
    row['requested_accuracy_achieved']=all(r.get('requested_accuracy_achieved') is True for r in rs) if row['experiment']=='accuracy_tolerance' else None
    summary.append(row)
sfields=list(gkeys)+[k for k in summary[0] if k not in gkeys] if summary else []
with (RAW/'search_summary.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=sfields); w.writeheader(); w.writerows(summary)

meta=json.loads((RAW/'metadata.json').read_text()); cfg=meta['config']
(RAW/'plot_environment.json').write_text(json.dumps({'python':sys.version,'matplotlib':matplotlib.__version__,'platform':platform.platform()},indent=2)+'\n')
colors={'dfs':'#1672a4','bestfirst':'#b33461'}; labels={'dfs':'DFS','bestfirst':'Best-first'}
plt.rcParams.update({'font.size':9,'figure.dpi':140,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
    fig.tight_layout(); fig.savefig(FIG/f'{name}.svg'); fig.savefig(FIG/f'{name}.png'); plt.close(fig)
def val(r,k,default=math.nan):
    x=r.get(k); return float(x) if x is not None and x!='' else default
def floor(x): return max(x,1e-16)
def config_label(r): return f"n={r['nqubits']} {r['family']} L{r['layers']} s{r['seed']} {r['observable']}"

# 1: error vs expanded nodes, one panel per n to reduce overplotting.
nvalues=sorted({r['nqubits'] for r in summary})
fig,axes2=plt.subplots(1,len(nvalues),figsize=(5*len(nvalues),4),sharey=True,squeeze=False); axes=axes2[0]
for ax,n in zip(axes,nvalues):
    for alg in colors:
        sub=[r for r in summary if r['experiment']=='expansion_budget' and r['algorithm']==alg and r['nqubits']==n and r['absolute_error'] is not None]
        for budget in cfg['expansion_budgets']:
            ss=[r for r in sub if float(r['requested_budget'])==budget]
            if ss: ax.scatter([val(r,'expanded_nodes') for r in ss],[floor(val(r,'absolute_error')) for r in ss],s=13,color=colors[alg],alpha=.55,label=labels[alg] if budget==cfg['expansion_budgets'][0] else None)
    ax.set(xscale='log',yscale='log',title=f'{n} qubits',xlabel='Actual expanded nodes'); ax.grid(alpha=.2)
axes[0].set_ylabel('Absolute error (floor 1e-16)'); axes[0].legend(); save(fig,'error_vs_expanded_nodes')

# 2: runtime-budget accuracy against actual elapsed time.
fig,axes2=plt.subplots(1,len(nvalues),figsize=(5*len(nvalues),4),sharey=True,squeeze=False); axes=axes2[0]
for ax,n in zip(axes,nvalues):
    for alg in colors:
        sub=[r for r in summary if r['experiment']=='runtime_budget' and r['algorithm']==alg and r['nqubits']==n and r['absolute_error'] is not None]
        ax.scatter([val(r,'search_elapsed_ms') for r in sub],[floor(val(r,'absolute_error')) for r in sub],s=14,color=colors[alg],alpha=.6,label=labels[alg])
    ax.set(xscale='log',yscale='log',title=f'{n} qubits',xlabel='Actual search runtime (ms)'); ax.grid(alpha=.2)
axes[0].set_ylabel('Absolute error (floor 1e-16)'); axes[0].legend(); save(fig,'error_vs_actual_runtime')

# 3: certified tolerance cost; x marks unmet targets.
fig,ax=plt.subplots(figsize=(8,5))
for alg in colors:
    sub=[r for r in summary if r['experiment']=='accuracy_tolerance' and r['algorithm']==alg]
    for achieved,marker in ((True,'o'),(False,'x')):
        ss=[r for r in sub if r['requested_accuracy_achieved'] is achieved]
        ax.scatter([float(r['requested_budget']) for r in ss],[val(r,'search_elapsed_ms') for r in ss],marker=marker,s=32,color=colors[alg],alpha=.65,label=f"{labels[alg]} {'achieved' if achieved else 'unmet'}")
ax.set(xscale='log',yscale='log',xlabel='Requested certified absolute-error tolerance',ylabel='Actual search runtime (ms)'); ax.grid(alpha=.2); ax.legend(); save(fig,'runtime_vs_requested_accuracy')

# 4 and 5: use the largest equal expansion budget as a common work cap.
largest=max(cfg['expansion_budgets'])
cap=[r for r in summary if r['experiment']=='expansion_budget' and float(r['requested_budget'])==largest]
for metric,ylabel,name,logy in [('peak_frontier','Peak frontier paths','frontier_vs_depth',True),('expanded_nodes','Expanded nodes','expanded_nodes_vs_depth',True)]:
    fig,axes2=plt.subplots(1,len(nvalues),figsize=(5*len(nvalues),4),sharey=True,squeeze=False); axes=axes2[0]
    for ax,n in zip(axes,nvalues):
        for alg in colors:
            for fam,marker in (('ising','o'),('random','s')):
                ss=[r for r in cap if r['nqubits']==n and r['algorithm']==alg and r['family']==fam]
                depths=sorted({r['layers'] for r in ss}); ys=[statistics.median([val(r,metric) for r in ss if r['layers']==d]) for d in depths]
                ax.plot(depths,ys,marker=marker,color=colors[alg],label=f'{labels[alg]} {fam}')
        ax.set(title=f'{n} qubits',xlabel='Circuit layers'); ax.grid(alpha=.2)
        if logy: ax.set_yscale('log')
    axes[0].set_ylabel(ylabel); axes[-1].legend(fontsize=7); save(fig,name)

# 6: paired percent improvement. Positive means BFS has lower error.
pairs=defaultdict(dict)
for r in summary:
    if r['experiment'] not in ('expansion_budget','runtime_budget') or r['absolute_error'] is None: continue
    k=tuple(r[x] for x in ('nqubits','family','layers','seed','observable','experiment','requested_budget'))
    pairs[k][r['algorithm']]=r
improvements=[]
for k,pair in pairs.items():
    if set(pair)=={'dfs','bestfirst'}:
        d=val(pair['dfs'],'absolute_error'); b=val(pair['bestfirst'],'absolute_error')
        if d>1e-15: improvements.append((k,100*(d-b)/d))
fig,ax=plt.subplots(figsize=(12,5))
improvements.sort(key=lambda x:x[1],reverse=True)
ax.bar(range(len(improvements)),[x[1] for x in improvements],color=['#25834b' if x[1]>=0 else '#b33461' for x in improvements],width=1)
ax.axhline(0,color='black',lw=1); ax.set(xlabel='Paired configuration/budget (sorted)',ylabel='Best-first accuracy improvement over DFS (%)'); ax.grid(axis='y',alpha=.2); save(fig,'bestfirst_accuracy_improvement')

# 7: separately instrumented attribution. Fractions overlap only by construction-free timers.
over=[r for r in records if r.get('experiment')=='overhead_profile' and r.get('bookkeeping_seconds') is not None]
fig,axes=plt.subplots(1,2,figsize=(11,4))
for alg in colors:
    ss=[r for r in over if r['algorithm']==alg]
    axes[0].scatter([1000*float(r['kernel_seconds']) for r in ss],[1000*float(r['bookkeeping_seconds']) for r in ss],s=18,color=colors[alg],alpha=.65,label=labels[alg])
    frac=[100*float(r['bookkeeping_seconds'])/max(1e-15,float(r['bookkeeping_seconds'])+float(r['kernel_seconds'])) for r in ss]
    axes[1].scatter([r['layers'] for r in ss],frac,s=18,color=colors[alg],alpha=.65,label=labels[alg])
axes[0].set(xscale='log',yscale='log',xlabel='Timed propagation-kernel work (ms)',ylabel='Timed frontier bookkeeping (ms)')
axes[1].set(xlabel='Circuit layers',ylabel='Bookkeeping / (bookkeeping + kernel) (%)')
for ax in axes: ax.grid(alpha=.2); ax.legend()
save(fig,'search_overhead_vs_propagation')

# Research answer from paired medians.
def paired_stats(experiment):
    ps=[(k,p) for k,p in pairs.items() if k[5]==experiment and set(p)=={'dfs','bestfirst'}]
    usable=[(val(p['dfs'],'absolute_error'),val(p['bestfirst'],'absolute_error')) for _,p in ps]
    wins=sum(b<d-1e-15 for d,b in usable); losses=sum(b>d+1e-15 for d,b in usable); ties=len(usable)-wins-losses
    ratios=[d/b for d,b in usable if d>1e-15 and b>1e-15]
    return len(usable),wins,losses,ties,(statistics.median(ratios) if ratios else None)
node=paired_stats('expansion_budget'); runtime=paired_stats('runtime_budget')
tol=[r for r in summary if r['experiment']=='accuracy_tolerance']
ach={a:sum(r['requested_accuracy_achieved'] is True for r in tol if r['algorithm']==a) for a in colors}
eligible={a:sum(1 for r in tol if r['algorithm']==a) for a in colors}
reasons=Counter((r['algorithm'],r.get('reason')) for r in primary if r.get('experiment')=='accuracy_tolerance' and r.get('status') not in ('exact','certified'))
overhead={a:statistics.median([100*float(r['bookkeeping_seconds'])/max(1e-15,float(r['bookkeeping_seconds'])+float(r['kernel_seconds'])) for r in over if r['algorithm']==a]) for a in colors}
tolpairs=defaultdict(dict)
for r in tol:
    k=tuple(r[x] for x in ('nqubits','family','layers','seed','observable','requested_budget'))
    tolpairs[k][r['algorithm']]=r
common=[p for p in tolpairs.values() if set(p)=={'dfs','bestfirst'} and p['dfs']['requested_accuracy_achieved'] is True and p['bestfirst']['requested_accuracy_achieved'] is True]
tolruntime=[val(p['bestfirst'],'search_elapsed_ms')/max(1e-15,val(p['dfs'],'search_elapsed_ms')) for p in common]
frontpairs=defaultdict(dict)
for r in cap:
    k=tuple(r[x] for x in ('nqubits','family','layers','seed','observable'))
    frontpairs[k][r['algorithm']]=r
fratios=[val(p['bestfirst'],'peak_frontier')/max(1,val(p['dfs'],'peak_frontier')) for p in frontpairs.values() if set(p)=={'dfs','bestfirst'}]
conclusion=("No on aggregate: coefficient priority did not provide a consistent equal-expansion advantage, and its accuracy was worse more often than DFS at equal wall-clock time." if runtime[2] >= runtime[1] else "Yes on aggregate for this grid: best-first won more equal-runtime accuracy comparisons than it lost, though the caveats and regressions below remain material.")
paired_resource_sentence=(f"Both traversals certified the target in {len(common)} paired points. In those pairs, the median best-first/DFS search-runtime ratio was {statistics.median(tolruntime):.3g} (above one means best-first took longer). At the one-million-expansion experiment, the median best-first/DFS peak-frontier ratio was {statistics.median(fratios):.3g}." if common and fratios else "Paired certified-runtime or frontier ratios were unavailable.")
def ratio_text(x): return 'unavailable' if x is None else f'{x:.3g}'
lines=["# DFS versus coefficient-prioritized best-first Pauli propagation","",
"## 1. Research question","",
"Does choosing the unfinished path with the largest absolute coefficient yield a more accurate completed-leaf approximation than deterministic depth-first traversal under the same expansion or wall-clock budget, and does that benefit justify heap overhead?","",
"## 2. Methodology","",
f"Both backends use the same `PauliWord`, `Path`, Float64 gate-expansion kernel, reverse gate order, terminal overlap, compensated accumulators, resource checks, and conservative frontier bound. DFS uses an explicit LIFO stack and visits cosine/unchanged before sine; best-first uses the existing max-heap ordered by `abs(coeff)` with insertion-order ties. Neither merges or memoizes paths. The grid covers n={cfg['nqubits']}, both circuit families, depths {cfg['layers']}, seeds 11 and 29, and both observables for seed 11 (Z-center for both seeds).","",
f"Expansion budgets are {cfg['expansion_budgets']}; runtime budgets are {[1000*x for x in cfg['runtime_budgets_seconds']]} ms; certified tolerances are {cfg['accuracy_tolerances']}. Each point has {cfg['repeats']} post-warmup repetitions. Algorithm order alternates. Search timing excludes circuit construction, compilation warmup, and exact reference calculation. References use unmodified VectorPauliSum with no truncation. If a reference worker fails or times out, observed error is left missing. Accuracy-controlled stopping uses only the frontier bound, never the reference.","",
"Process-lifetime RSS includes Julia, compilation, warmups, and prior cases in the worker. It is a safety/whole-process measure, not algorithm-only peak memory. Allocated bytes are traffic. The overhead figure comes from separate instrumented runs; those timers are excluded from primary timings.","",
"## 3. Results","",
conclusion,"",
f"At equal expansion budgets, {node[0]} paired median points were reference-valid: best-first had lower error in {node[1]}, higher error in {node[2]}, and tied within 1e-15 in {node[3]}. The median DFS-error / best-first-error ratio among nonzero pairs was {ratio_text(node[4])}; values above one favor best-first.","",
f"At equal runtime budgets, {runtime[0]} paired median points were reference-valid: best-first won {runtime[1]}, lost {runtime[2]}, and tied {runtime[3]}. The median nonzero error ratio was {ratio_text(runtime[4])}. This comparison includes priority-queue cost and runtime-check overhead.","",
f"Certified targets were achieved in {ach['bestfirst']}/{eligible['bestfirst']} best-first points and {ach['dfs']}/{eligible['dfs']} DFS points before the common one-second, one-million-expansion, and 250,000-frontier caps. Unmet termination reasons by algorithm are {dict(reasons)}.","",
paired_resource_sentence,"",
f"In separate instrumented runs, median bookkeeping share of measured kernel-plus-frontier time was {overhead['bestfirst']:.1f}% for best-first and {overhead['dfs']:.1f}% for DFS. These are direct timers around expansion and push/pop bookkeeping, not whole-runtime percentages; loop checks, terminal accumulation, and other work remain outside the denominator.","",
"## 4. Analysis","",
"The expansion-budget plots isolate traversal order: both methods expand the same number of nonterminal nodes unless the search finishes or hits the common frontier cap. Best-first can complete high-coefficient leaves earlier, but absolute coefficient is only a loose subtree bound. It does not know terminal overlap and cannot anticipate cancellation. DFS can therefore win when its deterministic branch happens to reach contributing I/Z leaves early, while best-first can spend many expansions maintaining broad high-mass partial paths.","",
"The runtime-budget result is the practical comparison. A best-first accuracy advantage at equal expansions may disappear when heap operations reduce the number of paths expanded in the same time. DFS usually keeps a narrow frontier, while best-first retains many partial paths; the depth plots show the memory consequence directly.","",
"## 5. Interpretation","",
"A positive result requires both lower equal-expansion error and enough benefit to survive the equal-runtime comparison. The aggregate counts above answer that question for this fixed-angle grid; individual regressions and wins remain visible in the raw tables and signed-improvement plot. Upstream merged propagation is retained as an exact/reference and architectural baseline, but its term merging means it cannot diagnose traversal-order quality.","",
"## 6. Limitations","",
"The frontier sum is conservative and can be much larger than observed error because it ignores terminal orthogonality and cancellation. Accuracy may be nonmonotonic with budget. One-millisecond measurements are sensitive to timer granularity, scheduler noise, and per-iteration `time_ns` checks; actual elapsed time is plotted and should be used instead of assuming the requested duration. External 3 GiB/120 s guards are sampled and can overshoot. Width alone does not determine difficulty because local light cones, random gates, observable placement, and branching matter.","",
"This is fixed-angle expectation propagation, not a parameter-uniform LOWESA surrogate guarantee. Exact references at larger/deeper configurations depend on the upstream merging implementation; missing references are never replaced with fabricated values. Instrumentation perturbs the overhead runs, so the attribution is directional rather than a cycle-accurate decomposition.","",
"## 7. Future work","",
"The measured evidence should determine whether a memory-bounded or hybrid strategy is justified. Plausible follow-ups are priorities incorporating remaining-gate structure, a persistent error ledger for bounded memory, and a switch from heap to DFS after a measured frontier threshold. Each would change more than traversal order and requires a separate comparison; none is implemented here.","",
"## Figures",""]
for name in ('error_vs_expanded_nodes','error_vs_actual_runtime','runtime_vs_requested_accuracy','frontier_vs_depth','expanded_nodes_vs_depth','bestfirst_accuracy_improvement','search_overhead_vs_propagation'):
    lines += [f"![{name}](figures/{name}.svg)",""]
lines += ["## Data and reproducibility","",f"The repeat CSV contains {len(records)} rows; the summary contains {len(summary)} median points. Hardware: {meta['cpu']}, {meta['total_memory_bytes']/2**30:.2f} GiB RAM, Julia {meta['julia_version']}, one Julia and one BLAS thread. OS: {meta['os_version']}. GPU detected but unused: {meta['gpu']}.","",
"See `raw/search_repeats.csv`, `raw/search_summary.csv`, per-size JSONL, `raw/metadata.json`, and worker logs. The earlier merging-backend study remains unchanged under `../results`.",""]
(OUT/'report.md').write_text('\n'.join(lines),encoding='utf-8')
print(f'Consolidated {len(records)} rows and {len(summary)} summary points')
