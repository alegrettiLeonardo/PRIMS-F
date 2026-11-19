"""Frozen-time 3-equation tangent experiment, not a new nonlinear NSCBC solver.

Includes the actual source Jacobian. Characteristic variant retains outgoing
increments and imposes the SAME mass-flux/pressure conditions as the current
operator, without interpreting R^-1 U as global invariants.
"""
from pathlib import Path
import json
import numpy as np
from scipy.sparse import csc_matrix,eye
from scipy.sparse.linalg import splu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_local import Audit,fd,scalar
from audit_paths import OUT
def main():
 a=Audit(OUT/'release/audit.so');u=a.state();q=a.point(u);a.par[17:19]=q['f'][:2];a.par[19]=q['f'][0]/q['v'][1];a.setup()
 A=q['a'];S=q['ds']/a.par[3]
 eig,R=np.linalg.eig(A);assert np.max(abs(eig.imag))==0;order=np.argsort(eig);eig=eig[order].real;R=R[:,order].real;L=np.linalg.inv(R)
 alpha,rl,rg=q['v'][:3]
 den=rl*alpha/a.par[13]**2+rg*(1-alpha)/a.par[12]**2
 pg=np.array([rg,rl,0])/den*a.par[2]
 D=R@np.diag(abs(eig))@L;fl=.5*(A+D);fr=.5*(A-D)
 J={};B={};details={}
 for side in [1,2]:
  current=fd(lambda x:a.boundary(x,side=side)['f'],u,1e-4)
  DB=fd(lambda x:a.boundary(x,side=side)['bc'],u,1e-4)
  incoming=eig>0 if side==1 else eig<0;outgoing=~incoming
  C=A[:2] if side==1 else pg[None,:]
  left=L[outgoing,:];T=np.vstack([C,left]);rhs=np.vstack([np.zeros((C.shape[0],3)),left])
  assert T.shape==(3,3)
  char=np.linalg.solve(T,rhs)
  J[side]={'current':current,'characteristic':A@char};B[side]={'current':DB,'characteristic':char}
  details[str(side)]={'eigenvalues_m_s':eig,'constraint_rank':int(np.linalg.matrix_rank(C@R[:,incoming])),
    'constraint_count':len(C),'incoming_count':int(sum(incoming)),
    'current_outgoing_increment_defect':left@(DB-np.eye(3)),
    'characteristic_outgoing_increment_defect':left@(char-np.eye(3)),
    'current_constraint_derivative':C@DB,'characteristic_constraint_derivative':C@char}
 n=80;dx=1.;dt=.002;duration=3.;positions=(np.arange(n)+.5)*dx
 output={'meaning':'frozen-time tangent about a nonstationary uniform state; actual source Jacobian included; not nonlinear/transient qualification',
         'parameters':a.par,'state_nd':u,'dt_s':dt,'dx_m':dx,'details':details,'runs':{}}
 fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,sharey=True)
 for row,direction in enumerate(['left','right']):
  wave=0 if direction=='left' else 2
  acoustic=R[:,wave]/(pg@R[:,wave]);initial=np.outer(np.exp(-((positions-40)/4)**2)*100,acoustic).ravel()
  for col,variant in enumerate(['current','characteristic']):
   M=np.kron(np.eye(n),S)
   def block(i,j,m):M[3*i:3*i+3,3*j:3*j+3]+=m
   for i in range(n-1):
    block(i,i,-fl/dx);block(i,i+1,-fr/dx);block(i+1,i,fl/dx);block(i+1,i+1,fr/dx)
   block(0,0,J[1][variant]/dx);block(n-1,n-1,-J[2][variant]/dx)
   lu=splu(eye(3*n,format='csc')-dt*csc_matrix(M));x=initial.copy();frames=[];times=[]
   for step in range(round(duration/dt)+1):
    if step%10==0:frames.append(x.reshape(n,3).copy());times.append(step*dt)
    x=lu.solve(x)
   data=np.array(frames);pressure=np.einsum('tnk,k->tn',data,pg)
   modal=np.einsum('jk,tnk->tnj',L,data)*np.array(pg@R)[None,None,:]
   np.savez_compressed(OUT/f'acoustic_{direction}_{variant}.npz',time_s=times,x_m=positions,delta_U=data,pressure_Pa=pressure,modal_pressure_Pa=modal)
   reflected=2 if direction=='left' else 0
   output['runs'][f'{direction}_{variant}']={'max_pressure_Pa':float(np.max(abs(pressure))),
    'reflected_modal_peak_Pa':float(np.max(abs(modal[np.array(times)>.9,:,reflected]))),
    'source_jacobian_included':True}
   axes[row,col].pcolormesh(positions,times,pressure,shading='auto',cmap='RdBu_r',vmin=-100,vmax=100)
   axes[row,col].set_title(direction+' / '+variant);axes[row,col].set_xlabel('s (m)');axes[row,col].set_ylabel('t (s)')
 fig.suptitle('Perturbação de pressão (Pa): operador tangente congelado; fontes incluídas')
 fig.tight_layout();fig.savefig(OUT/'acoustic_tangent.png',dpi=160);plt.close(fig)
 (OUT/'acoustic_tangent.json').write_text(json.dumps(scalar(output),indent=2)+'\n')
 print(json.dumps(output['runs'],indent=2))
if __name__=='__main__':main()
