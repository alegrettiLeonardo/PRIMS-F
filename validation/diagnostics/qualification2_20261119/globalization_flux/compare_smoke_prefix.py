from pathlib import Path
import json

root=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119')
new=root/'active_set_smokes_v1/D_SS1'
old=Path('/home/leonardoab/Documentos/PRIMS-F-numerics-results-20261119/v1_2400/runs/D_SS1')
records={}
for name in ('monitor.csv','fields.csv'):
    a=(old/name).read_text().splitlines()
    b=(new/name).read_text().splitlines()
    aa=[a[0]]+[row for row in a[1:] if float(row.split(',')[0])<78.5]
    bb=[b[0]]+[row for row in b[1:] if float(row.split(',')[0])<78.5]
    n=min(len(aa),len(bb))
    records[name]={'baseline_rows':len(aa),'smoke_rows':len(bb),'compared_rows':n,'equal_prefix':aa[:n]==bb[:n]}
(root/'offline_replay/default_prefix_check.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2))
