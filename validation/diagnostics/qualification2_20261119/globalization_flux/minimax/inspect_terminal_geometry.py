import argparse
import json
from offline_probe import Probe,InvalidState,ROOT,np

parser=argparse.ArgumentParser()
parser.add_argument('--case',default='D_SS3_active_terminal')
args=parser.parse_args()
p=Probe(args.case)
x=p.trial.copy()
r,j,n,_,_=p.evaluate(x)
s=1e-8*p.scale
a=j*s[None,:]/s[:,None]
b=r/s
index=17
coordinates=range(198,201)
gradient=np.zeros(306)
ulp=[]
for k in coordinates:
    h=1e-5*max(abs(x[k]),1e-5 if k%3==1 else 1e-2)
    values=[]
    for offset in (-2,-1,1,2):
        y=x.copy();y[k]+=offset*h
        values.append(p.margins(y)[index])
    gradient[k]=(values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
    item={'variable':k,'U':float(x[k]),'gradient':float(gradient[k]),'neighbors':[]}
    for direction in (-np.inf,np.inf):
        y=x.copy();y[k]=np.nextafter(x[k],direction)
        outcome={'delta_U':float(y[k]-x[k]),'margin_si':float(p.margins(y)[index])}
        try:
            outcome['residual_norm']=p.evaluate(y,False)[2]
            outcome['status']=0
        except InvalidState as error:
            outcome.update(error.details)
        item['neighbors'].append(outcome)
    ulp.append(item)
g=p.margins(x)[index]
c=gradient*s/1e-8
z0=np.linalg.solve(a,-b)
dual=np.linalg.solve(a.T,c)
delta=-g/1e-8-c@z0
ls_e=dual*delta/(dual@dual)
ls_z=z0+np.linalg.solve(a,ls_e)
inf_bound=abs(delta)/np.sum(np.abs(dual))
inf_e=np.sign(dual)*delta/np.sum(np.abs(dual))
inf_z=z0+np.linalg.solve(a,inf_e)
out={'case':args.case,'initial_norm':n,'endpoint_margin_si':float(g),'condition_scaled_J_2':float(np.linalg.cond(a)),
     'linear_ls_optimum_inf_norm':float(np.max(np.abs(ls_e))*1e-8),
     'linear_ls_optimum_l2_norm':float(np.linalg.norm(ls_e)*1e-8),
     'linear_minimum_possible_inf_norm_single_constraint':float(inf_bound*1e-8),
     'q_l1':float(np.sum(np.abs(dual))),'q_l2':float(np.linalg.norm(dual)),
     'state_ulp_probes':ulp,'proposals':[]}
for method,z in [('L2_KKT',ls_z),('LINF_LINEAR',inf_z)]:
    for exponent in (0,1,17):
        lam=2.**-exponent;y=x+lam*s*z
        item={'method':method,'lambda':lam,'endpoint_margin_si':float(p.margins(y)[index]),
              'max_delta_U_over_fixed_scale':float(np.max(np.abs(y-x)/p.scale))}
        try:
            item['residual_norm']=p.evaluate(y,False)[2]
            item['status']=0
        except InvalidState as error:
            item.update(error.details)
        out['proposals'].append(item)
np.savez(ROOT/(args.case+'_linear_geometry.npz'),x=x,r=r,j=j,scale=p.scale,gradient=gradient,
         unconstrained_step=s*z0,kkt_step=s*ls_z,linf_step=s*inf_z)
(ROOT/(args.case+'_linear_geometry.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
