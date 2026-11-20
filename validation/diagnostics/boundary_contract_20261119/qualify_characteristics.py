from pathlib import Path
import csv,json
import numpy as np
from audit_local import Audit,scalar
from audit_paths import OUT
def main():
 a=Audit(OUT/'variants_v2/standard/audit.so');p=a.par;p[0]=1.16;p[3]=p[4]/p[0];p[16]=206000/p[2];p[17]=545.38635/(p[1]*p[0]);p[18]=1.396675899947762/(p[1]*p[0]);p[19]=.547/p[0];a.models[:]=[1,3];a.setup()
 result={}
 for variant in ['standard','volume']:
  buckets={}
  with (OUT/'windows_v2'/variant/'boundary_states.csv').open() as f:
   for row in csv.DictReader(f):
    for side in [1,2]:
     for kind,pre in [('interior','Ui'),('boundary','UB')]:
      prefix=pre+('L' if side==1 else 'R');u=np.array([float(row[prefix+str(k)]) for k in [1,2,3]])
      theta=float(row['theta_in_rad' if side==1 else 'theta_out_rad']);q=a.point(u,theta,1 if side==1 else 3);assert q['status']==0
      vals,R=np.linalg.eig(q['a']);assert np.max(abs(vals.imag))==0;L=np.linalg.inv(R);incoming=vals>0 if side==1 else vals<0
      al,rl,rg=q['v'][:3];den=rl*al/p[13]**2+rg*(1-al)/p[12]**2;pg=np.array([rg,rl,0])/den
      C=q['a'][:2].copy() if side==1 else pg[None,:]
      Q=C.copy()
      if side==1:Q[0]=(q['a'][0]-q['f'][0]/rl*pg/p[12]**2)/rl
      key=f'{side}_{kind}';b=buckets.setdefault(key,{'frames':0,'incoming_counts':set(),'legacy_ranks':set(),'physical_ranks':set(),'source_projection_max_relative':0.,'min_eig_abs_m_s':float('inf')})
      b['frames']+=1;b['incoming_counts'].add(int(sum(incoming)))
      b['legacy_ranks'].add(int(np.linalg.matrix_rank((C/np.linalg.norm(C,axis=1)[:,None])@R[:,incoming])))
      b['physical_ranks'].add(int(np.linalg.matrix_rank((Q/np.linalg.norm(Q,axis=1)[:,None])@R[:,incoming])))
      b['source_projection_max_relative']=max(b['source_projection_max_relative'],np.max(abs(R@(L@q['s'])-q['s']))/max(1.,np.max(abs(q['s']))))
      b['min_eig_abs_m_s']=min(b['min_eig_abs_m_s'],float(np.min(abs(vals))*p[0]))
  for b in buckets.values():
   for k in ['incoming_counts','legacy_ranks','physical_ranks']:b[k]=sorted(b[k])
  result[variant]=buckets
 (OUT/'characteristic_ranks.json').write_text(json.dumps(scalar(result),indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
