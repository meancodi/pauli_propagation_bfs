import json
from pathlib import Path
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'results'/'figures'; FIG.mkdir(parents=True,exist_ok=True)
def read(name):
    with open(ROOT/'results'/name) as f: return [json.loads(x) for x in f if x.strip()]

s=read('synthetic.jsonl')
for pattern in ('local','random'):
    fig,ax=plt.subplots()
    for b in ('dictionary','dense','sparse'):
        rows=[x for x in s if x['pattern']==pattern and x['backend']==b]
        for occ in sorted({x['requested_occupancy'] for x in rows}):
            r=[x for x in rows if x['requested_occupancy']==occ]
            ax.plot([x['terms'] for x in r],[x['seconds'] for x in r],marker='o',label=f'{b}, occ={occ:g}')
    ax.set(xlabel='input Pauli terms',ylabel='gate time (s)',title=f'Fixed-region synthetic: {pattern} gates')
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(FIG/f'synthetic_{pattern}_runtime.png',dpi=150); plt.close(fig)

u=read('upstream_circuits.jsonl')
fig,ax=plt.subplots()
for b in ('PauliSum','VectorPauliSum','DenseLCB','SparseLCB'):
    r=[x for x in u if x['backend']==b]
    ax.plot([x['nqubits']+0.02*x['depth'] for x in r],[x['seconds'] for x in r],'.-',label=b)
ax.set(xlabel='width (depth offset)',ylabel='end-to-end time (s)',title='Upstream circuit comparison'); ax.legend(); fig.tight_layout(); fig.savefig(FIG/'upstream_runtime.png',dpi=150); plt.close(fig)

g=read('regrouping.jsonl') if (ROOT/'results'/'regrouping.jsonl').exists() else []
if g:
    fig,ax=plt.subplots()
    for b in ('DenseLCB','SparseLCB'):
        r=[x for x in g if x['backend']==b and x['terms']==10000]
        ax.plot([x['gates_per_region'] for x in r],[x['seconds'] for x in r],'.-',label=b)
    ax.set(xlabel='local gates per active-region segment',ylabel='8-segment time (s)',title='Regrouping break-even campaign'); ax.legend(); fig.tight_layout(); fig.savefig(FIG/'regrouping_break_even.png',dpi=150); plt.close(fig)
