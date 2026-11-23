import json
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results'/'validated_v2'; FIG=OUT/'figures'; FIG.mkdir(parents=True,exist_ok=True)
def read(name):
    p=OUT/name
    return [json.loads(x) for x in p.read_text().splitlines() if x.strip()] if p.exists() else []
rows=read('break_even.jsonl'); real=read('realistic_intermediate.jsonl')
complete=read('complete_circuits.jsonl')
done=[x for x in rows if x.get('status')=='complete']

def save(fig,name):
    fig.tight_layout(); fig.savefig(FIG/name,dpi=160); plt.close(fig)

for mode in sorted({x['mode'] for x in done}):
    for n in sorted({x['n'] for x in done}):
        for occ in sorted({x['requested_occupancy'] for x in done if x['n']==n and x['mode']==mode}):
            fig,ax=plt.subplots()
            for b,color in [('VectorPauliSum','tab:blue'),('PauliSum','tab:orange'),('DenseLCB','tab:green'),('SparseLCB','tab:red')]:
                r=[x for x in done if x['mode']==mode and x['n']==n and x['requested_occupancy']==occ and x['backend']==b and x['segments']==8]
                r=sorted(r,key=lambda x:x['gates_per_segment'])
                if r: ax.errorbar([x['gates_per_segment'] for x in r],[x['construction_inclusive']['median'] for x in r],
                    yerr=[x['construction_inclusive']['mad'] for x in r],marker='o',label=b,color=color)
            ax.set(xlabel='primitive gates per segment',ylabel='end-to-end runtime (s)',title=f'{mode}, n={n}, requested occupancy={occ:g}, 8 segments')
            ax.legend(fontsize=8); save(fig,f'runtime_{mode}_n{n}_occ{occ:g}.png')

# Paired DenseLCB speedup relative to VectorPauliSum, with zero at the center.
for mode in sorted({x['mode'] for x in done}):
    fig,ax=plt.subplots()
    keys=sorted({(x['n'],x['requested_occupancy']) for x in done if x['mode']==mode})
    for n,occ in keys:
        vals=[]
        for g in sorted({x['gates_per_segment'] for x in done}):
            v=[x for x in done if x['mode']==mode and x['n']==n and x['requested_occupancy']==occ and x['segments']==8 and x['gates_per_segment']==g]
            d=next((x for x in v if x['backend']=='DenseLCB'),None); u=next((x for x in v if x['backend']=='VectorPauliSum'),None)
            if d and u: vals.append((g,100*(u['construction_inclusive']['median']/d['construction_inclusive']['median']-1)))
        if vals: ax.plot([a for a,_ in vals],[b for _,b in vals],'.-',label=f'n={n}, occ={occ:g}')
    ax.axhline(0,color='black',lw=.8); ax.set(xlabel='primitive gates per segment',ylabel='DenseLCB improvement over VectorPauliSum (%)',title=f'Paired end-to-end speedup: {mode}')
    ax.legend(fontsize=7); save(fig,f'speedup_{mode}.png')

fig,ax=plt.subplots()
for b,c in [('DenseLCB','tab:green'),('SparseLCB','tab:red')]:
    r=[x for x in done if x['backend']==b and x['segments']==8 and x['gates_per_segment']==16]
    for n in sorted({x['n'] for x in r}):
        q=sorted([x for x in r if x['n']==n],key=lambda x:x['requested_occupancy'])
        if q: ax.plot([x['requested_occupancy'] for x in q],[x['allocated_bytes']/1e6 for x in q],'.-',label=f'{b}, n={n}',color=c,alpha=.7)
ax.set(xlabel='requested occupancy',ylabel='allocated bytes (MB)',title='Allocation versus occupancy, 8 segments, g=16'); ax.legend(fontsize=7); save(fig,'allocations_vs_occupancy.png')

fig,ax=plt.subplots()
for n in sorted({x['n'] for x in done}):
    r=[x for x in done if x['backend']=='DenseLCB' and x['n']==n and x['mode']=='alternating']
    r=sorted(r,key=lambda x:x['gates_per_segment'])
    if r: ax.plot([x['gates_per_segment'] for x in r],[x['regroup_seconds'] for x in r],'.-',label=f'n={n}')
ax.set(xlabel='primitive gates per segment',ylabel='regrouping time (s)',title='Regrouping component, alternating supports'); ax.legend(); save(fig,'regrouping_time.png')

if real:
    fig,ax=plt.subplots()
    labels=[]; v=[]; d=[]
    for x in real:
        labels.append(f"{x['family']}\nL{x['layers']}\n{x['seed']}"); v.append(x['vector_seconds']); d.append(x['dense_lcb_seconds'])
    idx=range(len(labels)); ax.bar([i-.18 for i in idx],v,.36,label='VectorPauliSum',color='tab:blue'); ax.bar([i+.18 for i in idx],d,.36,label='DenseLCB',color='tab:green')
    ax.set_xticks(list(idx),labels,rotation=45,ha='right'); ax.set_ylabel('continuation runtime (s)'); ax.set_title('Captured intermediate operators'); ax.legend(); save(fig,'realistic_intermediate_runtime.png')

for family in sorted({x['family'] for x in complete}):
    fig,ax=plt.subplots()
    for n in sorted({x['nqubits'] for x in complete if x['family']==family}):
        r=[x for x in complete if x['family']==family and x['nqubits']==n]
        r=sorted(r,key=lambda x:x['layers'])
        ax.plot([x['layers'] for x in r],[x['vector_seconds'] for x in r],'.-',label=f'VectorPauliSum n={n}')
        ax.plot([x['layers'] for x in r],[x['dense_lcb_seconds'] for x in r],'.--',label=f'DenseLCB n={n}')
    ax.set(xlabel='circuit layers',ylabel='runtime (s)',title=f'Complete circuit comparison: {family}'); ax.legend(fontsize=8); save(fig,f'complete_{family}_runtime.png')
