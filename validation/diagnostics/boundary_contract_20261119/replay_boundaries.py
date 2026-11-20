from pathlib import Path
import ctypes as ct, json, sys
import numpy as np
from scipy.linalg import solve
from audit_local import Audit,arr,ptr,I,D,fd,scalar
from audit_paths import ROOT
from audit_paths import OUT
def evaluate(lib,x,jac=0):
 r=arr(306);j=arr((306,306));norm=D();fin=arr(3);fout=arr(3);status=I();idx=I();stage=I()
 lib.probe_eval(ptr(x),I(jac),ptr(r),ptr(j),ct.byref(norm),ptr(fin),ptr(fout),ct.byref(status),ct.byref(idx),ct.byref(stage))
 return {'status':status.value,'index':idx.value,'stage':stage.value,'norm':norm.value,'r':r,'j':j,'fin':fin,'fout':fout}
def main():
 result={}
 arrays=OUT/'replay_arrays';arrays.mkdir(exist_ok=True)
 for variant in ['standard','volume','partition']:
  a=Audit(OUT/'variants'/variant/'audit.so');lib=a.lib;result[variant]={}
  for case in ['D_SS1','D_SS3','D_SS1_active_terminal','D_SS3_active_terminal']:
   folder=ROOT/'validation/diagnostics/qualification2_20261119/globalization_flux/fixtures'/case
   encoded=str(folder).encode();st=I();lib.probe_init(ct.c_char_p(encoded),I(len(encoded)),I(1),ct.byref(st));assert st.value==0
   values={}
   for which in [0,1]:
    x=arr(306);sc=arr(306);meta=arr(6);lib.probe_state(I(which),ptr(x),ptr(sc),ptr(meta))
    q=evaluate(lib,x,1);endpoints={}
    np.savez_compressed(arrays/(variant+'_'+case+'_'+str(which)+'.npz'),state_nd=x,residual_nd=q['r'],jacobian_nd=q['j'],scales_nd=sc,metadata=meta)
    snapshot=np.loadtxt(folder/'last_rejected_state.dat',skiprows=2)
    for side,k in [(1,0),(2,101)]:
     u=x[3*k:3*k+3];theta=snapshot[k,1];b=a.boundary(u,theta,side,1);p=a.point(b['bc'],theta,1 if side==1 else 3,1) if b['status']==0 else None
     endpoints[str(side)]={'boundary':b,'point':p}
     if p is not None and p['status']==0:
      ev=np.linalg.eigvals(p['a']);endpoints[str(side)]['eigenvalues_SI']=ev*meta[3]
      endpoints[str(side)]['incoming']=int(np.sum(ev.real>0 if side==1 else ev.real<0))
    values[str(which)]={'time_s':meta[0],'dt_s':meta[1],'operator':{k:v for k,v in q.items() if k not in ['r','j']},'endpoints':endpoints}
    if q['status']==0:
     if which==1:
      direction=solve(q['j'],-q['r']);proposals=[]
      for ls in range(18):
       lam=2.0**(-ls);candidate=x+lam*direction;ev=evaluate(lib,candidate)
       ends=[]
       for side,k in [(1,0),(2,101)]:
        bd=a.boundary(candidate[3*k:3*k+3],snapshot[k,1],side,1)
        ends.append({'side':side,'status':bd['status'],'flux':bd['f'],'p_face_nd':bd['p']})
       interface=arr(48);ist=I();lib.probe_interface(ptr(candidate),ptr(interface),ct.byref(ist))
       gas_trace=None
       if ist.value==0:
        rho_g_star=interface[3]*meta[5]/343**2/meta[4]
        gas_trace=interface[1]/rho_g_star*meta[3]
       proposals.append({'lambda':lam,'status':ev['status'],'index':ev['index'],'stage':ev['stage'],
                        'norm':ev['norm'],'boundary':ends,'conservative_candidate_nd':candidate,
                        'junction_probe_status':ist.value,'junction_data_nd':interface,'jg_star_m_s':gas_trace})
      values[str(which)]['newton_ray']=proposals
     # Full residual directional derivative with boundary reconstruction, same frozen old dissipation.
     # Direction supported only in endpoint variables; compare complete 306 residual components.
     v=np.zeros(306);v[:3]=sc[:3]*np.array([.37,-.23,.11]);v[-3:]=sc[-3:]*np.array([-.31,.19,.43])
     predicted=q['j']@v;curve=[]
     for h in [1e-3,1e-4,1e-5,1e-6,1e-7]:
      plus=evaluate(lib,x+h*v);minus=evaluate(lib,x-h*v)
      item={'h':h,'statuses':[plus['status'],minus['status']]}
      if plus['status']==minus['status']==0:
       obs=(plus['r']-minus['r'])/(2*h);item['scaled_error']=float(np.max(abs(obs-predicted)/sc))
      curve.append(item)
     values[str(which)]['full_boundary_direction_fd']=curve
   result[variant][case]=values
 (OUT/'replays.json').write_text(json.dumps(scalar(result),indent=2)+'\n')
 for v,cases in result.items():
  print(v,{c:{q:x['operator']['status'] for q,x in states.items()} for c,states in cases.items()},flush=True)
if __name__=='__main__':main()
