import json, pathlib, shutil

R=pathlib.Path(__file__).parents[1]/'results'/'hierarchical_v2'
src=R/'benchmark.jsonl'; backup=R/'benchmark.interleaved_backup.jsonl'; clean=R/'benchmark.repaired.jsonl'
if not backup.exists(): shutil.copy2(src,backup)
records={}; invalid=0
for line in src.open(encoding='utf-8'):
    try:
        x=json.loads(line)
        if isinstance(x,dict) and 'config_id' in x: records[x['config_id']]=x
    except json.JSONDecodeError: invalid+=1
with clean.open('w',encoding='utf-8',newline='\n') as f:
    for x in records.values(): f.write(json.dumps(x,separators=(',',':'))+'\n')
shutil.copy2(clean,src)
(R/'repair_manifest.json').write_text(json.dumps({'source':str(backup),'repaired':str(src),'valid_records':len(records),'invalid_lines_discarded':invalid},indent=2)+'\n')
print(f'repaired {len(records)} records; discarded {invalid} invalid/interleaved lines')
