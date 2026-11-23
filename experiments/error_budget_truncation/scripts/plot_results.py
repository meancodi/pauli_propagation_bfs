import json, pathlib
import matplotlib.pyplot as plt

root=pathlib.Path(__file__).resolve().parents[1]
rows=[json.loads(x) for x in (root/'results/raw/benchmark.jsonl').read_text().splitlines()]
budget=[r for r in rows if r['algorithm']=='budgeted']
if not budget: raise SystemExit('no benchmark data')
out=root/'results/figures'; out.mkdir(parents=True,exist_ok=True)
def grouped(field, ylabel, filename):
    fig, ax=plt.subplots(figsize=(8,5))
    for policy in sorted({r['policy'] for r in budget}):
        rr=[r for r in budget if r['policy']==policy and r['nqubits']==5]
        by={}
        for r in rr: by.setdefault(r['epsilon_total'],[]).append(r[field])
        xs=sorted(by); ys=[sum(v)/len(v) for x in xs for v in [by[x]]]
        ax.plot(xs,ys,'o-',label=policy)
    ax.set_xscale('symlog',linthresh=1e-4); ax.set_xlabel('global certified error budget'); ax.set_ylabel(ylabel); ax.legend(); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(out/filename,dpi=160); plt.close(fig)
grouped('elapsed_ms','runtime (ms)','runtime_vs_budget.png')
grouped('peak_terms','peak Pauli terms','peak_terms_vs_budget.png')
grouped('observed_error','absolute expectation error','actual_error_vs_budget.png')
