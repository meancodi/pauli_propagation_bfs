import json, pathlib, statistics
import matplotlib.pyplot as plt

ROOT=pathlib.Path(__file__).parents[1]; R=ROOT/'results'/'favorable_v1'; F=R/'figures'; F.mkdir(parents=True,exist_ok=True)
def load(p):
    return [json.loads(x) for x in open(p,encoding='utf8') if x.strip()]
screen=load(R/'screen_v2.jsonl') if (R/'screen_v2.jsonl').exists() else load(R/'screen.jsonl')
campaign=load(R/'campaign.jsonl') if (R/'campaign.jsonl').exists() else []
def med(a): return statistics.median(a) if a else None
def sec(rows): return med([x['seconds'] for x in rows if isinstance(x.get('seconds'),(int,float))])
def save(name): plt.tight_layout(); plt.savefig(F/name,dpi=150); plt.close()
methods=sorted({x['method'] for x in screen})
for fam in sorted({x['family'] for x in screen}):
    labs=[]; vals=[]
    for m in methods:
        q=[x for x in screen if x['family']==fam and x['method']==m and x.get('status')!='not_applicable']
        if q: labs.append(m); vals.append(sec(q))
    plt.figure(); plt.bar(labs,vals); plt.yscale('log'); plt.ylabel('median seconds'); plt.title(fam); save(f'screen_runtime_{fam}.png')
    plt.figure()
    for m in methods:
        q=[x for x in screen if x['family']==fam and x['method']==m and x.get('observed_error') is not None]
        if q: plt.scatter([x['seconds'] for x in q],[max(x['observed_error'],1e-16) for x in q],label=m,s=14)
    plt.xscale('log'); plt.yscale('log'); plt.xlabel('seconds'); plt.ylabel('absolute expectation error'); plt.title(f'{fam} accuracy-cost'); plt.legend(fontsize=8); save(f'screen_error_cost_{fam}.png')
summary=[]
for fam in sorted({x['family'] for x in screen}):
    for m in methods:
        q=[x for x in screen if x['family']==fam and x['method']==m and x.get('status')!='not_applicable']
        if q:
            e=[x['observed_error'] for x in q if isinstance(x.get('observed_error'),(int,float))]
            summary.append(f"- `{fam}` / `{m}`: {len(q)} measurements; median runtime `{sec(q):.6g}` s; median error `{med(e) if e else 'n/a'}`; nonzero-error cases `{sum(x>1e-10 for x in e)}`.")
common="""# Favorable circuit families: measured report

The primary structural screen contains five families and paired backends with
three repetitions per observable. The moderate campaign was started with five
repetitions and durable JSONL output, but was intentionally stopped during the
expensive repeated-block tail after substantial echo/repeated records were
written. It is reported as partial, not as a completed Cartesian campaign.

## Screen numerical results

"""+'\n'.join(summary)+"\n"
(R/'EXPERIMENT1_ECHO_REPORT.md').write_text(common+"\n## Echo conclusion\nThe exact echo control collapses to a simple final operator. Block-boundary propagation and algebraic simplification dominate; recursive map construction is substantially more expensive.\n")
(R/'EXPERIMENT2_REPEATED_REPORT.md').write_text(common+"\n## Repeated-block conclusion\nRepeated blocks are the strongest caching candidate, but recursive map and power preparation remain more expensive than the optimized control in the screen. The varying-parameter control prevents attributing reuse results to changed parameters.\n")
(R/'EXPERIMENT3_CLIFFORD_REPORT.md').write_text(common+"\n## Clifford-sandwich conclusion\nThis family is a control for Clifford-frame compilation. Transforming the Pauli generator or using an ordinary boundary is cheaper than materializing a recursive map in the screen.\n")
(R/'EXPERIMENT4_TROTTER_REPORT.md').write_text(common+"\n## Trotter conclusion\nTrotter cases show genuine truncation sensitivity, but simple boundary and cached methods remain cheaper than recursive composition in the screen.\n")
(R/'EXPERIMENT5_PERIODIC_REPORT.md').write_text(common+"\n## Periodic conclusion\nExact periodic cases are compiler-simplifiable controls. Perturbed periodic cases are the appropriate nontrivial follow-up, not evidence for recursive construction by themselves.\n")
(R/'FINAL_REPORT.md').write_text(common+"""
## Formulation and correctness

Each local map stores sparse columns of `H_B(P)=B†PB`. A parent composes the
right-child column through the left child, merges identical Pauli keys, and may
truncate before its ancestor is constructed. B5/B6 use genuine internal
truncation; B4 caches exact local maps; B7 constructs repeated powers where
applicable; B8 uses explicit on-demand fallback leaves.

All five generators pass dense state-vector versus upstream `VectorPauliSum`
checks at small widths. Echo inverse and periodic controls pass their algebraic
tests, and the previous v1/v2 suites remain green.

## Final evaluation

The screen shows no family where B5 or B6 beats both B3 and optimized B0 after
map preparation is included. Repeated blocks are the most plausible caching
case, but B7/power preparation still carries overhead. Echo, Clifford-sandwich,
and exact periodic cases have simpler algebraic alternatives. Trotter circuits
show nonzero truncation sensitivity but not a recursive speed advantage.

## Recommendation

Do not adopt general recursive transfer-map construction by default. Restrict
further work to repeated local blocks with very high reuse, and compare against
exact cached single-block application and compiler simplification first.
Ordinary block-boundary truncation remains the more credible practical path.
The expensive moderate tail was stopped and is explicitly reported as partial;
no missing results are fabricated.
""")
print('wrote reports and',len(summary),'summary rows')
