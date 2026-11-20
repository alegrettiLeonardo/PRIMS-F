"""Attach hashes and terminal solver metrics without modifying physical data."""
from pathlib import Path
import hashlib,json,subprocess
from audit_paths import ROOT,OUT
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 index=[]
 for group in ['windows','windows_v2','controls']:
  for run in sorted((OUT/group).iterdir()):
   if not run.is_dir() or not (run/'run_status.txt').is_file():continue
   variant=run.name if group.startswith('windows') else run.name.rsplit('_',1)[-1]
   build=OUT/('variants_v2' if group=='windows_v2' else 'variants')/variant
   inp=run/'input.dat' if group=='controls' else OUT/'windows/OSC_frame3000_input.dat'
   status=dict(line.split('=',1) for line in (run/'run_status.txt').read_text().splitlines())
   record={'run':str(run),'operator_base':'2c3e3af8c99d490ffc1dcbf9357da6c0d59d633d','variant':variant,
     'build_manifest':str(build/'manifest.json'),'build_manifest_sha256':sha(build/'manifest.json'),
     'binary_sha256':sha(build/'window'),'input':str(inp),'input_sha256':sha(inp),
     'terminal_status':{k:v.strip() for k,v in status.items()},
     'compiler':'GNU Fortran 13.3.0','environment':{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},
     'outputs':{p.name:sha(p) for p in run.iterdir() if p.is_file() and p.name!='audit_manifest.json'},
     'limits':['Experimental variant identified by generated source hashes; production files unchanged.',
               'OSC initial condition is reconstructed SI data, not an exact historic restart.'] if group.startswith('windows') else ['Declared synthetic boundary control, not an experimental slugging case.']}
   (run/'audit_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
   index.append({'run':str(run.relative_to(OUT)),'manifest_sha256':sha(run/'audit_manifest.json'),
                 'status':record['terminal_status']['status'],'reached_time_s':record['terminal_status']['reached_time_s']})
 (OUT/'run_manifest_index.json').write_text(json.dumps(index,indent=2)+'\n');print(len(index),'run manifests')
if __name__=='__main__':main()
