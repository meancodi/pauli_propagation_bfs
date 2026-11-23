from pathlib import Path
from collections import defaultdict
import csv,json,math,itertools
root=Path(__file__).resolve().parents[1]; out=root/'results_dfs_vs_bestfirst'; raw=out/'raw'
rows=list(csv.DictReader((raw/'search_repeats.csv').open(newline='',encoding='utf-8')))
assert len(rows)==len({r['run_id'] for r in rows})
for r in rows:
    if r.get('absolute_error'):
        e=float(r['absolute_error']); b=float(r['frontier_bound'])
        assert math.isfinite(e) and e<=b+1e-10, r['run_id']
    if r.get('status')=='exact':
        assert float(r['frontier_bound'])==0 and int(r['frontier_nodes'])==0
    if r.get('algorithm') in ('dfs','bestfirst') and r.get('peak_frontier'):
        assert int(r['peak_frontier'])<=250000
# Equal-expansion pairs must use identical counts unless one completed or hit frontier first.
groups=defaultdict(dict)
for r in rows:
    if r.get('experiment')=='expansion_budget':
        k=tuple(r[x] for x in ('nqubits','family','layers','seed','observable','requested_budget','repeat'))
        groups[k][r['algorithm']]=r
for k,p in groups.items():
    if set(p)=={'dfs','bestfirst'}:
        d,b=p['dfs'],p['bestfirst']; budget=int(float(d['requested_budget']))
        for r in (d,b):
            assert int(r['expanded_nodes'])<=budget
            if r['status']=='resource_limited' and r['reason']=='expansions': assert int(r['expanded_nodes'])==budget
cfg=json.loads((raw/'metadata.json').read_text())['config']
for name in ('error_vs_expanded_nodes','error_vs_actual_runtime','runtime_vs_requested_accuracy','frontier_vs_depth','expanded_nodes_vs_depth','bestfirst_accuracy_improvement','search_overhead_vs_propagation'):
    assert (out/'figures'/f'{name}.svg').is_file()
print(f'Validated {len(rows)} DFS/BFS records and {len(groups)} equal-expansion pairs')
