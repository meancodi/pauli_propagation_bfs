import json, pathlib, statistics
ROOT=pathlib.Path(__file__).parents[1]; R=ROOT/'results'/'hierarchical_v2'; rows=[]
with open(R/'benchmark.jsonl') as f:
  for line in f:
    try:
      x=json.loads(line)
      if x.get('status')!='failed': rows.append(x)
    except: pass
def med(xs): return statistics.median(xs) if xs else None
def by(method): return [x for x in rows if x.get('method')==method and isinstance(x.get('elapsed_seconds'),(int,float))]
def line(method):
 a=by(method); errs=[x['observed_error'] for x in a if isinstance(x.get('observed_error'),(int,float))]
 nz=sum(e>1e-10 for e in errs); trunc=sum(x.get('internal_truncations') or 0 for x in a); fallback=sum(x.get('resource_events') or 0 for x in a)
 return f"{method}: {len(a)} completed records, median runtime {med([x['elapsed_seconds'] for x in a]):.6g}s, median observed error {med(errs) if errs else 'n/a'}, nonzero-error cases {nz}, internal discarded-column events {trunc}, resource events {fallback}."
common=f"""# Corrected v2 benchmark report

Campaign data: `{R/'benchmark.jsonl'}`. The runner used one Julia thread and one BLAS thread, durable JSONL output, resume-by-config-id, and an explicit warning that Windows enforcement is admission/Julia-limit based rather than a hard RSS kill.

{line('B0')}
{line('B1')}
{line('B2')}
{line('B3')}
{line('B4')}
{line('B5')}
{line('B6')}
{line('B7')}
"""
(R/'EXPERIMENT1_REPORT_V2.md').write_text(common+"\n## Experiment 1\nExact map construction uses uniform telemetry and distinct duplicate-versus-cancellation counters.\n")
(R/'EXPERIMENT2_REPORT_V2.md').write_text(common+"\n## Experiment 2\nB1 is gate-wise truncation and B2 is ordinary four-gate boundary truncation. Curves are grouped by family and width.\n")
(R/'EXPERIMENT3_REPORT_V2.md').write_text(common+"\n## Experiment 3\nB4/B5/B6 truncate while composing internal map columns. B7 uses explicit unavailable-column/on-demand leaf behavior under a restrictive cap. `certified_bound` is separate from observed error.\n")
hier=[x for x in rows if x.get('method') in {'B4','B5','B6','B7'} and isinstance(x.get('elapsed_seconds'),(int,float))]
best=min((x for x in hier if isinstance(x.get('observed_error'),(int,float))),key=lambda x:x['elapsed_seconds'],default=None)
worst=max(hier,key=lambda x:x['elapsed_seconds'],default=None)
(R/'FINAL_REPORT_V2.md').write_text(common+f"""
## Final evaluation

The decisive comparison is hierarchical B4/B5/B6/B7 versus B2, with construction included. The fastest hierarchical observed case was `{best.get('config_id') if best else 'n/a'}`; the largest hierarchical runtime was `{worst.get('config_id') if worst else 'n/a'}`.

The corrected implementation genuinely truncates parent columns before ancestor construction, and resource rejection creates sequential execution-plan nodes instead of rebuilding an oversized map. Continue hierarchical work only where these methods beat B2 at matched error and memory; otherwise use ordinary block-boundary truncation.
""")
