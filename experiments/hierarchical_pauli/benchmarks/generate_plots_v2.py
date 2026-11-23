import json, pathlib, statistics
import matplotlib.pyplot as plt

ROOT=pathlib.Path(__file__).parents[1]; R=ROOT/'results'/'hierarchical_v2'; F=R/'figures'; F.mkdir(parents=True,exist_ok=True)
rows=[]
with open(R/'benchmark.jsonl') as fh:
    for line in fh:
        try:
            x=json.loads(line)
            if x.get('status') not in ('failed','skipped'): rows.append(x)
        except json.JSONDecodeError: pass
def val(x,k): return x.get(k) if isinstance(x.get(k),(int,float)) else None
def save(name): plt.tight_layout(); plt.savefig(F/name,dpi=150); plt.close()
methods=sorted({x['method'] for x in rows}); colors={m:('green' if m in {'B4','B5','B6','B7'} else 'tab:blue') for m in methods}
for family in sorted({x['family'] for x in rows}):
  for n in sorted({x['nqubits'] for x in rows if x['family']==family}):
    plt.figure()
    for m in methods:
      a=[x for x in rows if x['family']==family and x['nqubits']==n and x['method']==m and val(x,'elapsed_seconds') is not None]
      if a: plt.plot([x['layers'] for x in a],[x['elapsed_seconds'] for x in a],'.',label=m,color=colors[m],alpha=.7)
    plt.yscale('log'); plt.xlabel('layers'); plt.ylabel('end-to-end seconds'); plt.title(f'{family}, n={n}'); plt.legend(ncol=2,fontsize=8); save(f'runtime_{family}_n{n}.png')
for n in sorted({x['nqubits'] for x in rows}):
  plt.figure()
  for m in methods:
    a=[x for x in rows if x['nqubits']==n and x['method']==m and val(x,'observed_error') is not None and val(x,'elapsed_seconds') is not None]
    if a: plt.scatter([x['elapsed_seconds'] for x in a],[max(x['observed_error'],1e-16) for x in a],s=8,label=m,color=colors[m])
  plt.xscale('log'); plt.yscale('log'); plt.xlabel('seconds'); plt.ylabel('absolute expectation error'); plt.title(f'error-cost, n={n}'); plt.legend(ncol=2,fontsize=8); save(f'error_vs_runtime_n{n}.png')
for family in sorted({x['family'] for x in rows}):
  for n in sorted({x['nqubits'] for x in rows if x['family']==family}):
    vals=[]
    for m in methods:
      if m=='B2': continue
      a=[x for x in rows if x['family']==family and x['nqubits']==n and x['method']==m and val(x,'elapsed_seconds') is not None]
      b=[x for x in rows if x['family']==family and x['nqubits']==n and x['method']=='B2' and val(x,'elapsed_seconds') is not None]
      if a and b: vals.append((m,100*(statistics.median(x['elapsed_seconds'] for x in b)/statistics.median(x['elapsed_seconds'] for x in a)-1)))
    if vals:
      plt.figure(); names=[x[0] for x in vals]; ys=[x[1] for x in vals]; plt.axhline(0,color='black'); plt.bar(names,ys,color=['green' if y>=0 else 'red' for y in ys]); plt.ylabel('improvement over B2 (%)'); plt.title(f'{family}, n={n}'); save(f'improvement_over_B2_{family}_n{n}.png')
