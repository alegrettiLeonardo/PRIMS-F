"""Small independent stationary solve of the actual finite-volume residual.

Full finite differences include re-evaluation of the stationary dissipation.
No replacement residual is used on rejected points. Target 1e-11, max20 Newton,
18 halvings, set before execution. Not the production transient algorithm.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import solve
from audit_local import Audit,fd,scalar
from audit_paths import OUT
def main():
 results={}
 for variant in ['standard','volume']:
  a=Audit(OUT/'variants_v2'/variant/'audit.so');n=4;theta=np.full(n,-.035);width=np.full(n,.1)
  u=np.repeat(a.boundary(a.state(),theta=-.035)['bc'][None,:],n,axis=0);trace=[]
  def evaluate(x):
   q=x.reshape(n,3);return a.assemble(q,old=q,theta=theta,width=width,npip=2,dt=1.)
  x=u.ravel();scale=np.tile([1.,.001,1.],n)
  for it in range(20):
   q=evaluate(x)
   if q['status']!=0:raise RuntimeError(('inadmissible base',q['status']))
   r=q['r'].ravel(order='F');norm=np.max(abs(r)/scale)
   trace.append({'iteration':it,'norm':norm})
   if norm<=1e-11:break
   def function(z):
    qz=evaluate(z)
    if qz['status']!=0:raise RuntimeError(('FD outside model domain',qz['status']))
    return qz['r'].ravel(order='F')
   J=fd(function,x,1e-4);dx=solve(J,-r);accepted=False
   for ls in range(18):
    lam=2.**(-ls);trial=evaluate(x+lam*dx)
    if trial['status']==0 and np.max(abs(trial['r'].ravel(order='F'))/scale)<norm:
     x=x+lam*dx;trace[-1]['lambda']=lam;accepted=True;break
   if not accepted:raise RuntimeError('stationary globalization exhausted')
  else:raise RuntimeError('stationary Newton exhausted')
  u=x.reshape(n,3);inlet=a.boundary(u[0],theta=-.035);outlet=a.boundary(u[-1],theta=-.035,side=2)
  rho=1+(inlet['p']-1)/a.par[12]**2
  results[variant]={'parameters':a.par,'theta':theta,'width':width,'state_nd':u,'trace':trace,
   'norm':norm,'inlet':inlet,'outlet':outlet,'Ql_relative_error':inlet['f'][0]/rho/a.par[19]-1,
   'p_out_error_Pa':(outlet['p']-a.par[16])*a.par[2],
   'net_phase_flux_nd':q['fin'][:2]-q['fout'][:2]}
 (OUT/'stationary_nonzero.json').write_text(json.dumps(scalar(results),indent=2)+'\n')
 print({k:{key:v[key] for key in ['norm','Ql_relative_error','p_out_error_Pa','net_phase_flux_nd']} for k,v in results.items()})
if __name__=='__main__':main()
