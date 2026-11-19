"""Bounded paired windows. Reconstructed printed fields are NOT exact old restarts."""
from pathlib import Path
import csv, hashlib, json, os, subprocess, time, math, sys
import numpy as np
from audit_local import Audit,scalar
from audit_paths import ROOT
from audit_paths import OUT
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 base=Path('/home/leonardoab/Documentos/PRIMS-F-numerics-results-20261119/v1_osc3600/runs/D_OSC')
 fields=base/'fields.csv';rows=[]
 with fields.open() as f:
  for row in csv.DictReader(f):
   t=float(row['time_s'])
   if t>=3000 and not rows:target=t
   if t>=3000 and t==target:rows.append({k:float(v) for k,v in row.items()})
   if rows and t>target:break
 assert len(rows)==102,len(rows)
 a=Audit(OUT/'release/audit.so');p=a.par.copy();p[0]=1.16;p[3]=p[4]/p[0];p[16]=206000/p[2]
 p[17]=545.38635/(p[1]*p[0]);p[18]=1.396675899947762/(p[1]*p[0]);p[19]=.547/p[0]
 states=[]
 for i,r in enumerate(rows):
  al=r['alpha'];rl=r['rhol_kg_m3']/p[1];rg=r['rhog_kg_m3']/p[1]
  u=[(1-al)*rl,al*rg,((1-al)*rl*r['ul_m_s']+al*rg*r['ug_m_s'])/p[0]]
  width=(p[5]/p[4] if i<51 else 1)/51
  states.append((r['theta_rad'],width,*u))
 states=[states[0],states[0],*states,states[-1],states[-1]]
 folder=OUT/'windows';folder.mkdir(exist_ok=True)
 lines=['106 2 51 1 3 120.0 1e-8 1e-8 25 F 3000',' '.join(format(x,'.17e') for x in p),
        ' '.join(format(x/p[3],'.17e') for x in [.001,.125,1e-7])]
 lines += [str(i+1)+' '+' '.join(format(x,'.17e') for x in r) for i,r in enumerate(states)]
 inp=folder/'OSC_frame3000_input.dat';inp.write_text('\n'.join(lines)+'\n')
 summary={'frame_time_s':target,'initial_state':'U reconstructed from binary64-printed SI primitives; common to all variants; not exact historic restart',
          'source_fields':str(fields),'source_fields_sha256':sha(fields),'input_sha256':sha(inp),'runs':[]}
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 for variant in ['plain','standard','volume']:
  dest=folder/variant;dest.mkdir(exist_ok=False)
  cmd=[str(OUT/'variants'/variant/'window'),str(inp),str(dest)]
  start=time.perf_counter()
  with (dest/'stdout.log').open('w') as f:
   try:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=180);code=r.returncode
   except subprocess.TimeoutExpired:code='TIMEOUT'
  item={'variant':variant,'command':cmd,'returncode':code,'wall_s':time.perf_counter()-start,
        'binary_sha256':sha(Path(cmd[0])),'input_sha256':sha(inp)}
  status=dest/'run_status.txt'
  if status.exists():item['run_status']=status.read_text()
  item['output_hashes']={f.name:sha(f) for f in dest.iterdir() if f.is_file()}
  summary['runs'].append(item);(folder/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
  print(variant,code,item['wall_s'],item.get('run_status','NO_FINAL_STATUS')[:220],flush=True)
 for name in ['monitor.csv','fields.csv','conservation.csv','newton_trials.csv']:
  x=folder/'plain'/name;y=folder/'standard'/name
  if x.exists() and y.exists():summary.setdefault('observer_identity',{})[name]=sha(x)==sha(y)
 (folder/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
