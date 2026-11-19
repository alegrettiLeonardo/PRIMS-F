#!/usr/bin/env python3
"""Archive completed/failed comparison evidence; no executables or bibliography."""
from pathlib import Path
import argparse,hashlib,json,subprocess,zipfile
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--repo',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=a.repo.resolve()
assert not a.output.exists(),a.output
folders=sorted((root/'validation/results').glob('comparison*'))
assert folders
manifests=[json.loads((q/'manifest.json').read_text()) for q in folders]
assert all(m['status'] in ['COMPLETED','FAILED','TIMEOUT'] for m in manifests)
files={}
allowed={'.txt','.json','.csv','.dat','.nml','.log','.md','.py','.f90','.png','.pdf'}
for folder in folders+[root/'validation/comparisons']:
  for q in sorted(folder.rglob('*')):
    if q.is_file() and q.suffix.lower() in allowed:files[str(q.relative_to(root))]=q
for m in manifests:
  deck=Path(m['command'][-1])
  if not deck.is_absolute():deck=root/deck
  assert hashlib.sha256(deck.read_bytes()).hexdigest()==m['namelist_sha256']
  files['effective_namelists/'+deck.name]=deck
for filename in ['run_validation.py','run_closure_comparison.py','analyze_validation.py','analyze_osc_envelopes.py','analyze_phase_sequence.py','validation_evidence.py','osc_run_selection.py','compare_validation.py']:
  files['scripts/'+filename]=root/'scripts'/filename
freeze=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
hashes={key:hashlib.sha256(q.read_bytes()).hexdigest() for key,q in files.items()}
manifest={'git_head':freeze,'runs':len(manifests),'statuses':{s:sum(m['status']==s for m in manifests) for s in ['COMPLETED','FAILED','TIMEOUT']},
 'contents_sha256':hashes,'scope':'raw numerical data, effective namelists, analyses, provenance and analysis scripts; no executables or bibliographic PDFs'}
a.output.parent.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(a.output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for key,q in files.items():z.write(q,key)
  z.writestr('CONTENTS_SHA256.json',json.dumps(manifest,indent=2)+'\n')
  z.writestr('README.txt','PRIMS-F numerical evidence archive\nSource commit: '+freeze+'\nUse CONTENTS_SHA256.json to verify individual files.\nFailed runs remain failed; pressure oscillations alone do not validate a regime.\n')
result={'path':str(a.output),'sha256':hashlib.sha256(a.output.read_bytes()).hexdigest(),'bytes':a.output.stat().st_size,
        'git_head':freeze,'files':len(files)+2,'runs':len(manifests),'statuses':manifest['statuses']}
a.output.with_suffix('.manifest.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
