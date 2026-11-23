import json, pathlib
import matplotlib.pyplot as plt

ROOT=pathlib.Path(__file__).parents[1]; OUT=ROOT/'results'/'figures'; OUT.mkdir(parents=True,exist_ok=True)
def load(name):
    with open(ROOT/'results'/name) as f: return [json.loads(x) for x in f if x.strip()]
e1,e2,e3=load('experiment1.jsonl'),load('experiment2.jsonl'),load('experiment3.jsonl')
def save(name): plt.tight_layout(); plt.savefig(OUT/name,dpi=150); plt.close()
for n in sorted({r['nqubits'] for r in e1}):
    plt.figure()
    for strategy in sorted({r['strategy'] for r in e1}):
        rs=[r for r in e1 if r['nqubits']==n and r['strategy']==strategy]
        plt.scatter([r['block_length'] for r in rs],[r['entries'] for r in rs],label=strategy)
    plt.xlabel('block length'); plt.ylabel('nonzero map entries'); plt.title(f'n={n}'); plt.legend(); save(f'experiment1_entries_n{n}.png')
    plt.figure()
    for strategy in sorted({r['strategy'] for r in e1}):
        rs=[r for r in e1 if r['nqubits']==n and r['strategy']==strategy]
        plt.scatter([r['block_length'] for r in rs],[r['seconds'] for r in rs],label=strategy)
    plt.yscale('log'); plt.xlabel('block length'); plt.ylabel('construction seconds'); plt.title(f'n={n}'); plt.legend(); save(f'experiment1_runtime_n{n}.png')
for n in sorted({r['nqubits'] for r in e2}):
    for field, ylabel, suffix in [('seconds','runtime (ms)','runtime'),('peak_terms','peak terms','peak_terms')]:
        plt.figure()
        for threshold in sorted({r['threshold'] for r in e2}):
            rs=[r for r in e2 if r['nqubits']==n and r['threshold']==threshold]
            plt.plot([r['block_size'] for r in rs],[r[field]*(1000 if field=='seconds' else 1) for r in rs],marker='o',label=str(threshold))
        plt.xlabel('block size'); plt.ylabel(ylabel); plt.title(f'n={n}'); plt.legend(title='threshold'); save(f'experiment2_{suffix}_n{n}.png')
for n in sorted({r['nqubits'] for r in e3}):
    rs=[r for r in e3 if r['nqubits']==n]; plt.bar([r['strategy'] for r in rs],[r['seconds']*1000 for r in rs]); plt.ylabel('end-to-end runtime (ms)'); plt.title(f'n={n}'); save(f'experiment3_runtime_n{n}.png')

# Corrected v2 plots are generated separately so the original smoke figures
# remain reproducible: python generate_plots_v2.py
