from pathlib import Path
import subprocess,time,hashlib,json,os,sys
from audit_local import Audit,scalar
import numpy as np
import csv
from audit_paths import ROOT
from audit_paths import OUT
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 here=Path(__file__).parent
 for debug in [False,True]:
  cmd=['python3',str(here/'build_variants.py'),'--root',str(ROOT),'--out',str(OUT/('variants_v2_debug' if debug else 'variants_v2')),'--variants','standard','volume','partition']
  if debug:cmd+=['--debug']
  subprocess.run(cmd,check=True)
 destroot=OUT/'windows_v2';destroot.mkdir(exist_ok=True);inp=OUT/'windows/OSC_frame3000_input.dat';runs=[]
 for v in ['standard','volume']:
  dest=destroot/v;dest.mkdir(exist_ok=False);cmd=[str(OUT/'variants_v2'/v/'window'),str(inp),str(dest)];start=time.perf_counter()
  with (dest/'stdout.log').open('w') as f:
   try:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=180,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));code=r.returncode
   except subprocess.TimeoutExpired:code='TIMEOUT'
  record={'variant':v,'command':cmd,'input_sha256':sha(inp),'binary_sha256':sha(cmd[0]),'returncode':code,'wall_s':time.perf_counter()-start,'outputs':{f.name:sha(f) for f in dest.iterdir() if f.is_file()}}
  record['observer_v1_identity']={n:sha(dest/n)==sha(OUT/'windows'/v/n) for n in ['monitor.csv','fields.csv','conservation.csv','boundary_contract.csv']}
  runs.append(record);(destroot/'manifest.json').write_text(json.dumps(runs,indent=2)+'\n');print('v2',v,code,record['observer_v1_identity'],flush=True)
 # Recompute point spectra from EXACT accepted U, never from interpolated primitives.
 a=Audit(OUT/'variants_v2/standard/audit.so');p=a.par;p[0]=1.16;p[3]=p[4]/p[0];p[16]=206000/p[2];p[17]=545.38635/(p[1]*p[0]);p[18]=1.396675899947762/(p[1]*p[0]);p[19]=.547/p[0];a.models[:]=[1,3];a.setup()
 stats={}
 for v in ['standard','volume']:
  samples=[];counts={}
  with (destroot/v/'boundary_states.csv').open() as f:
   for r in csv.DictReader(f):
    for side in [1,2]:
     for kind,prefix in [('interior','Ui'),('boundary','UB')]:
      prefix+='L' if side==1 else 'R';u=np.array([float(r[prefix+str(k)]) for k in [1,2,3]]);angle=float(r['theta_in_rad' if side==1 else 'theta_out_rad']);q=a.point(u,angle,1 if side==1 else 3)
      item={'step':int(r['step']),'time_s':float(r['time_s']),'side':side,'kind':kind,'status':q['status']}
      if q['status']==0:
       eig=np.linalg.eigvals(q['a']);item['real']=bool(np.max(abs(eig.imag))==0)
       item['incoming']=int(sum(eig.real>0 if side==1 else eig.real<0));item['zero']=int(sum(eig.real==0));item['eigenvalues_m_s']=eig.real*p[0];item['ul_m_s']=q['v'][4]*p[0];item['ug_m_s']=q['v'][5]*p[0]
      signature=(side,kind,item.get('status'),item.get('incoming'),item.get('zero'),item.get('real'));key=str(signature);counts[key]=counts.get(key,0)+1;samples.append(item)
  (destroot/v/'characteristics.json').write_text(json.dumps(scalar(samples),indent=2)+'\n');stats[v]=counts
 (destroot/'characteristic_counts.json').write_text(json.dumps(stats,indent=2)+'\n');print(stats)
if __name__=='__main__':main()
