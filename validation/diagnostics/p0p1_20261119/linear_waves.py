#!/usr/bin/env python3
"""Isolated frozen-linear wave verification; never a production solver selector.

Periodic unit interval. Each characteristic wave is advanced separately.
MUSCL/SSPRK2 varies only space; BE/SDIRK2 varies only time at fixed FV1.
The actual nonlinear primitive recovery qualifies every reconstructed trace.
"""
import argparse,ctypes as C,json,hashlib,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--fixture',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
 lib=C.CDLL(str(a.library.resolve()));ptr=C.POINTER(C.c_double);iptr=C.POINTER(C.c_int)
 base=np.zeros(3);A=np.zeros((3,3),order='F');v=np.zeros(6);theta=C.c_double();scales=np.zeros(3);st=C.c_int();s=str(a.fixture.resolve()).encode()
 lib.wave_setup.argtypes=[C.c_char_p,C.c_int,ptr,ptr,ptr,C.POINTER(C.c_double),ptr,iptr]
 lib.wave_setup(s,len(s),base.ctypes.data_as(ptr),A.ctypes.data_as(ptr),v.ctypes.data_as(ptr),C.byref(theta),scales.ctypes.data_as(ptr),C.byref(st));assert st.value==0
 lam,R=np.linalg.eig(A);assert np.max(np.abs(lam.imag))==0
 order=np.argsort(lam.real);lam=lam.real[order];R=R.real[:,order];Ri=np.linalg.inv(R)
 assert np.linalg.cond(R)<1e10 and np.max(np.abs(A@R-R*lam))<1e-9
 lib.wave_admissible.argtypes=[ptr,C.c_int,iptr]
 def admissible(u):
  u=np.asfortranarray(u,dtype=np.float64);status=np.zeros(len(u),dtype=np.int32);lib.wave_admissible(u.ctypes.data_as(ptr),len(u),status.ctypes.data_as(iptr));return status==0
 amplitudes=np.array([1e-3*min(base[k]/abs(R[k,i]) for k in [0,1] if R[k,i]!=0) for i in range(3)])
 linear=[];spatial=[];temporal=[];waveforms=[];limiter_uses=0
 for mode,speed in enumerate(lam):
  assert speed!=0
  for scheme in ['ROE','RUSANOV']:
   diffusion=abs(speed) if scheme=='ROE' else max(abs(lam))
   for n in [64,128,256]:
    h=1/n;k=2*np.pi/n;T=1/abs(speed);steps=4*n;dt=T/steps
    symbol=-1j*speed*np.sin(k)/h-diffusion*(1-np.cos(k))/h
    gain=(1-dt*symbol)**(-steps)
    q=np.sin(2*np.pi*(np.arange(n)+.5)/n)*np.sinc(1/n)
    flux=.5*(speed*(q+np.roll(q,-1))-diffusion*(np.roll(q,-1)-q));rhs=-(flux-np.roll(flux,1))/h
    measured=np.fft.rfft(rhs)[1]/np.fft.rfft(q)[1]
    assert abs(measured-symbol)<1e-10*max(1,abs(symbol)) and abs(rhs.sum())<1e-10
    linear.append(dict(mode=mode,speed_m_s=speed*scales[0],scheme=scheme,n=n,actual_dt_nd=dt,unit='prototype length/time in nondimensional characteristic units',one_period_gain=abs(gain),phase_error_rad=float(np.angle(gain)),artificial_diffusion_nd=diffusion))
  # Fixed semidiscrete FV1 symbol: errors below exclude spatial error.
  n=256;h=1/n;k=2*np.pi/n;T=1/abs(speed);symbol=-1j*speed*np.sin(k)/h-abs(speed)*(1-np.cos(k))/h;exact=np.exp(symbol*T)
  for method in ['BE','SDIRK2']:
   errors=[]
   for steps in [16,32,64,128]:
    dt=T/steps;z=dt*symbol;gamma=1-1/np.sqrt(2)
    gain=1/(1-z) if method=='BE' else (1+(1-2*gamma)*z)/(1-gamma*z)**2
    error=abs(gain**steps-exact);errors.append(error)
    temporal.append(dict(mode=mode,method=method,steps=steps,error_semidiscrete=error,actual_dt=dt))
   observed=math.log(errors[-2]/errors[-1],2);assert observed>(.8 if method=='BE' else 1.8)
  for n in [64,128,256]:
   x=(np.arange(n)+.5)/n;initial=np.sin(2*np.pi*x)*np.sinc(1/n);T=1/abs(speed);steps=5*n;dt=T/steps
   for space in ['FV1','MUSCL_MC']:
    q=initial.copy();eps=amplitudes[mode]
    def rhs(q):
     nonlocal limiter_uses
     slope=np.zeros(n)
     if space=='MUSCL_MC':
      dl=q-np.roll(q,1);dr=np.roll(q,-1)-q;dc=.5*(dl+dr)
      same=(np.sign(dl)==np.sign(dr))&(np.sign(dl)==np.sign(dc))
      slope=np.where(same,np.sign(dl)*np.minimum(np.minimum(2*np.abs(dl),2*np.abs(dr)),np.abs(dc)),0)
      assert np.all(admissible(base+eps*q[:,None]*R[:,mode]))
      for repeat in range(32):
       ok=admissible(base+eps*(q-.5*slope)[:,None]*R[:,mode])&admissible(base+eps*(q+.5*slope)[:,None]*R[:,mode])
       if ok.all():break
       # Mean-preserving trace limitation, never clipping conserved cell masses.
       slope[~ok]*=.5;limiter_uses+=int((~ok).sum())
      else:raise RuntimeError('No admissible reconstruction')
     ql=q+.5*slope;qr=np.roll(q-.5*slope,-1)
     flux=.5*(speed*(ql+qr)-abs(speed)*(qr-ql))
     return -(flux-np.roll(flux,1))*n
    for step in range(steps):
     stage=q+dt*rhs(q);q=.5*(q+stage+dt*rhs(stage))
    assert abs(q.sum()-initial.sum())<1e-9
    assert np.all(admissible(base+eps*q[:,None]*R[:,mode]))
    assert q.min()>=initial.min()-1e-12 and q.max()<=initial.max()+1e-12
    error=np.mean(abs(q-initial));spatial.append(dict(mode=mode,space=space,n=n,steps=steps,actual_dt=dt,l1_error=error,amplitude_gain=abs(np.fft.rfft(q)[1]/np.fft.rfft(initial)[1]),mass_sum_error=abs(q.sum()-initial.sum())))
    if mode==1 and n==128:waveforms.append((space,x.copy(),q.copy()))
 # Independent constraint activation: steep slopes have inadmissible endpoints.
 for mode in range(3):
  slope=R[:,mode]*min(base[:2]/np.maximum(abs(R[:2,mode]),1e-100))*20
  assert not admissible(np.array([base-slope/2,base+slope/2])).all()
  scale=1.
  for k in range(32):
   if admissible(np.array([base-scale*slope/2,base+scale*slope/2])).all():break
   scale*=.5
  assert 0<scale<1 and np.max(abs((base-scale*slope/2+base+scale*slope/2)/2-base))<1e-12
  limiter_uses+=k
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for scheme in ['ROE','RUSANOV']:
  d=[x for x in linear if x['mode']==1 and x['scheme']==scheme];axes[0].plot([x['n']for x in d],[x['one_period_gain']for x in d],'o-',label=scheme)
 axes[0].set(xlabel='volumes',ylabel='material wave gain after one period');axes[0].legend()
 for name,x,q in waveforms:axes[1].plot(x,q,label=name)
 axes[1].plot(x,np.sin(2*np.pi*x)*np.sinc(1/len(x)),'k--',label='exact cell average');axes[1].legend();fig.tight_layout();fig.savefig(a.output/'linear_waves.png',dpi=160);plt.close(fig)
 result=dict(scope='isolated linearization at one admissible physical state, independent periodic characteristic waves; no nonlinear severe-slugging qualification',library_sha256=hashlib.sha256(a.library.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),base_U=base.tolist(),base_primitives=v.tolist(),theta=theta.value,A=A.tolist(),eigenvalues_nd=lam.tolist(),scales=scales.tolist(),linear=linear,spatial=spatial,temporal=temporal,admissibility_slope_halvings=limiter_uses)
 (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS: physical spectrum, conservative FV1/Rusanov, MUSCL/TVD traces, BE/SDIRK2 temporal order; isolated only')
if __name__=='__main__':main()
