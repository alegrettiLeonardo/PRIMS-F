"""Independent boundary/operator probes, using the real Fortran routines via C ABI."""
from pathlib import Path
import argparse, ctypes as ct, json, math, os, subprocess, sys
import numpy as np

D=ct.c_double; I=ct.c_int
def ptr(x):return x.ctypes.data_as(ct.c_void_p)
def arr(shape):return np.zeros(shape,dtype=np.float64,order='F')

class Audit:
 def __init__(self,lib):
  self.lib=ct.CDLL(str(lib));self.par=np.array([1.,997.05,101325.,12.546253993197219,12.546253993197219,57.4,.05254,4.6e-5,9.81,.00089,.0000185,.0728,1498*math.sqrt(997.05/101325),343*math.sqrt(997.05/101325),1.,1.,2.1,.5,.001,.5]);self.models=np.array([2,2],dtype=np.int32);self.setup()
 def setup(self):self.lib.ba_setup(ptr(self.par),ptr(self.models))
 def boundary(self,u,theta=0.,side=1,replay=0):
  u=np.ascontiguousarray(u,dtype=float);f=arr(3);bc=arr(3);p=D();st=I()
  self.lib.ba_boundary(ptr(u),D(theta),I(side),I(replay),ptr(f),ptr(bc),ct.byref(p),ct.byref(st))
  return dict(status=st.value,f=f,bc=bc,p=p.value)
 def point(self,u,theta=0.,model=2,replay=0):
  u=np.ascontiguousarray(u,dtype=float);v=arr(6);f=arr(3);s=arr(3);a=arr((3,3));ds=arr((3,3));st=I()
  self.lib.ba_point(ptr(u),D(theta),I(model),I(replay),ptr(v),ptr(f),ptr(s),ptr(a),ptr(ds),ct.byref(st))
  return dict(status=st.value,v=v,f=f,s=s,a=a,ds=ds)
 def assemble(self,u,old=None,theta=None,width=None,dt=.001,npip=2,jac=0):
  u=np.asfortranarray(u,dtype=float);old=u.copy(order='F') if old is None else np.asfortranarray(old)
  n=len(u);theta=np.zeros(n) if theta is None else np.asarray(theta,dtype=float);width=np.ones(n) if width is None else np.asarray(width,dtype=float)
  r=arr((3,n));j=arr((3*n,3*n));fin=arr(3);fout=arr(3);base=arr((2,2));st=I();idx=I()
  self.lib.ba_assemble(I(n),I(npip),ptr(u),ptr(old),ptr(theta),ptr(width),D(dt),I(jac),ptr(r),ptr(j),ptr(fin),ptr(fout),ptr(base),ct.byref(st),ct.byref(idx))
  return dict(status=st.value,index=idx.value,r=r,j=j,fin=fin,fout=fout,base=base)
 def state(self,p=2.1,alpha=.25,mom=.6):
  rl=1+(p-1)/self.par[12]**2;rg=p/self.par[13]**2
  return np.array([(1-alpha)*rl,alpha*rg,mom])

def scalar(x):
 if isinstance(x,np.ndarray):return [scalar(t) for t in x.tolist()]
 if isinstance(x,(list,tuple)):return [scalar(t) for t in x]
 if isinstance(x,dict):return {k:scalar(v) for k,v in x.items()}
 if isinstance(x,np.generic):return scalar(x.item())
 if isinstance(x,float) and not math.isfinite(x):return None
 return x

def fd(fun,x,rel,order=4):
 ans=[]
 for k in range(len(x)):
  h=rel*max(abs(x[k]),1e-5 if k%3==1 else 1e-2)
  if k%3<2:h=min(h,.2*x[k])
  e=np.zeros(len(x));e[k]=h
  if order==2:q=(fun(x+e)-fun(x-e))/(2*h)
  else:q=(-fun(x+2*e)+8*fun(x+e)-8*fun(x-e)+fun(x-2*e))/(12*h)
  ans.append(q)
 return np.array(ans).T

def local(a):
 result={};u=a.state();b=a.boundary(u);o=a.boundary(u,side=2)
 result['algebra']={'inlet':b,'outlet':o,'mass_scale_kg_s':a.par[1]*a.par[0]*math.pi*a.par[6]**2/4,
  'liquid_Q_relative_error':b['f'][0]/(1+(b['p']-1)/a.par[12]**2)/a.par[19]-1,
  'outlet_pressure_error_Pa':(o['p']-a.par[16])*a.par[2]}
 # Independent EOS and conservation-flux identity; no reuse of boundary construction.
 for side,data in [(1,b),(2,o)]:
  q=a.point(data['bc']);assert data['status']==q['status']==0
  assert np.max(abs(data['f']-q['f'])/np.maximum(abs(data['f']),1))<1e-10
  if side==2:assert abs((data['p']-a.par[16])*a.par[2])<1e-7
 assert np.array_equal(b['f'][:2],a.par[17:19])
 # Independently perturb conservative variables, reconstructing the entire boundary each time.
 records=[]
 for side in [1,2]:
  for rel in [1e-3,3e-4,1e-4,3e-5,1e-5,3e-6,1e-6,3e-7,1e-7]:
   j2=fd(lambda z:a.boundary(z,side=side)['f'],u,rel,2)
   j4=fd(lambda z:a.boundary(z,side=side)['f'],u,rel,4)
   records.append({'side':side,'relative_h':rel,'j2':j2,'j4':j4,'scaled_error':np.max(abs(j2-j4)/np.maximum(1,abs(j4)))})
 result['boundary_fd']=records
 old=u[None,:];base=a.assemble(old,npip=1,jac=1);assert base['status']==0
 records=[]
 for rel in [1e-3,1e-4,1e-5,1e-6,1e-7]:
  jf=fd(lambda x:a.assemble(x.reshape(1,3),old=old,npip=1)['r'].ravel(order='F'),u,rel)
  records.append({'relative_h':rel,'max_scaled_error':np.max(abs(jf-base['j'])/np.maximum(1,abs(jf)))})
 result['full_jacobian']=records
 # Identical physical mesh, geometry, closure and U; only bookkeeping split changes.
 n=6;states=np.array([a.state(p=2.1+.001*i,alpha=.25+.01*i,mom=.6+.003*i) for i in range(n)])
 part=[]
 for th in [0.,-.035]:
  left=a.assemble(states,theta=np.full(n,th),npip=2);right=a.assemble(states,theta=np.full(n,th),npip=4)
  part.append({'theta':th,'status':[left['status'],right['status']],
   'residual_max_abs':np.max(abs(left['r']-right['r'])),'residuals':[left['r'],right['r']],
   'phase_base_mismatch':[np.max(abs(x['base'][:,0]-x['base'][:,1])) for x in [left,right]]})
 result['partition']=part
 # Counts at both sides and imposed/interior states, including reversed and pure flows.
 chars=[]
 for model,angle in [(2,0.),(1,-.035),(3,.9)]:
  a.models[:]=model;a.setup()
  for mom in [.6,-.6,-5.,0.]:
   for alpha in [.25,1e-10,0.,1.]:
    state=a.state(alpha=alpha,mom=mom)
    for side in [1,2]:
     bc=a.boundary(state,theta=angle,side=side)
     for kind,x in [('interior',state),('boundary',bc['bc'])]:
      item={'model':model,'angle_rad':angle,'momentum_nd':mom,'alpha_input':alpha,'side':side,'kind':kind,'boundary_status':bc['status']}
      if kind=='boundary' and bc['status']!=0:chars.append(item);continue
      q=a.point(x,theta=angle,model=model);item['point_status']=q['status'];item['primitive']=q['v']
      if q['status']==0 and np.all(np.isfinite(q['a'])):
       eig=np.linalg.eigvals(q['a']);item['eigenvalues_m_s']=eig*a.par[0]
       item['real_spectrum']=bool(np.max(abs(eig.imag))<1e-9)
       if item['real_spectrum']:item['incoming']=int(np.sum(eig.real>0 if side==1 else eig.real<0));item['zero']=int(np.sum(eig.real==0))
      chars.append(item)
 result['characteristics']=chars
 a.models[:]=2;a.setup()
 # Stagnant, filling/draining, gas-absent and gas-near-absent boundary operators.
 result['limits']=[]
 for ml,mg,alpha,mom in [(0,0,0.,0.),(.5,0,0.,.5),(0,0,1.,0.),(0,.001,1.,.001),(.5,0,1e-10,.5),(0,0,.25,-.1)]:
  a.par[17:19]=ml,mg;a.setup();state=a.state(alpha=alpha,mom=mom)
  result['limits'].append({'ml_nd':ml,'mg_nd':mg,'alpha':alpha,'momentum_nd':mom,'inlet':a.boundary(state),'outlet':a.boundary(state,side=2)})
 return result

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--library',type=Path,required=True);ap.add_argument('--output',type=Path);ap.add_argument('--phase-child',action='store_true');args=ap.parse_args()
 a=Audit(args.library)
 if args.phase_child:
  ct.CDLL('libm.so.6').feenableexcept(1|4|8)
  a.par[18]=0;a.setup();u=a.state(alpha=0.,mom=.5)[None,:]
  print('Calling one-cell physical boundary Jacobian at exactly absent gas',flush=True)
  print(scalar(a.assemble(u,npip=1,jac=1)),flush=True);return
 result=local(a);args.output.write_text(json.dumps(scalar(result),indent=2,allow_nan=False)+'\n')
 print(json.dumps(scalar({k:v for k,v in result.items() if k in ['algebra','full_jacobian','partition']}),indent=2))
if __name__=='__main__':main()
