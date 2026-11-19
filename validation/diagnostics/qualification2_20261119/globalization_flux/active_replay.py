"""Constrained linear least-squares Newton replay, without changing its residual."""
import argparse
import json
from offline_probe import Probe, InvalidState, ROOT, np

parser = argparse.ArgumentParser()
parser.add_argument('--case',default='D_SS1')
parser.add_argument('--fd-factor',type=float,default=1.0)
args = parser.parse_args()
p = Probe(args.case)
x0 = p.trial.copy()
xscale = 1e-8*p.scale
rscale = 1e-8*p.scale
margin_index = 0 if args.case=='D_SS1' else 17
margin_scale = 1e-8
indices = list(range(3*(53-3),3*(55-3))) if args.case=='D_SS1' else list(range(3*(69-3),3*(70-3)))


def constraint_gradient(x):
    grad = np.zeros(306)
    details = []
    for k in indices:
        h = args.fd_factor*1e-5*max(abs(x[k]),1e-5 if k%3==1 else 1e-2)
        values=[]
        for offset in (-2,-1,1,2):
            y=x.copy()
            y[k]+=offset*h
            values.append(p.margins(y)[margin_index])
        derivative=(values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
        grad[k]=derivative*xscale[k]/margin_scale
        details.append({'variable':k,'h':h,'values':values,'derivative':derivative})
    return grad,details


report={'case':args.case,'fd_factor':args.fd_factor,'method':'linear equality active-set Gauss-Newton with exact-residual feasible line search',
        'max_iterations':30,'max_backtracks':30,'constraint_index':margin_index,
        'residual_scaling':'1e-8 * original fixed Uold component scales',
        'variable_scaling':'1e-8 * original fixed Uold component scales',
        'margin_scaling_si':margin_scale,'iterations':[],'invalid_proposals':[]}
x=x0.copy()
for iteration in range(30):
    r,j,n,fi,fo=p.evaluate(x)
    g=p.margins(x)[margin_index]
    a=j*xscale[None,:]/rscale[:,None]
    b=r/rscale
    grad,details=constraint_gradient(x)
    unconstrained=np.linalg.solve(a,-b)
    # Minimum ||A dz+b|| subject to c dz=-g. Solve in residual space,
    # avoiding normal equations and preserving all conservative variables.
    dual=np.linalg.solve(a.T,grad)
    predicted=g/margin_scale+grad@unconstrained
    if predicted>=0:
        dz=unconstrained
        multiplier=0.
    else:
        correction=(-g/margin_scale-grad@unconstrained)/(dual@dual)
        dz=unconstrained+np.linalg.solve(a,dual*correction)
        multiplier=correction
    item={'iteration':iteration,'norm':n,'objective':.5*float(b@b),'margin_si':g,
          'predicted_unconstrained_margin_scaled':float(predicted),
          'linear_constrained_inf_norm':float(np.max(np.abs(a@dz+b))),
          'multiplier':float(multiplier),'gradient':details if iteration==0 else None}
    report['iterations'].append(item)
    if n<=1e-8:
        report['termination']='ORIGINAL_RESIDUAL_CONVERGED'
        break
    accepted=False
    for backtrack in range(30):
        lam=2.**-backtrack
        y=x+lam*xscale*dz
        try:
            gy=p.margins(y)[margin_index]
            if gy<0:
                report['invalid_proposals'].append({'iteration':iteration,'backtrack':backtrack,'lambda':lam,'constraint_si':float(gy),'reason':'ADMISSIBILITY'})
                continue
            rr,_,nn,_,_=p.evaluate(y,False)
        except InvalidState as error:
            report['invalid_proposals'].append({'iteration':iteration,'backtrack':backtrack,'lambda':lam,**error.details})
            continue
        objective=.5*float((rr/rscale)@(rr/rscale))
        if objective<item['objective']:
            item['accepted_lambda']=lam
            item['accepted_norm']=nn
            item['accepted_margin_si']=float(gy)
            item['accepted_objective']=objective
            x=y
            accepted=True
            break
    if not accepted:
        report['termination']='NO_ADMISSIBLE_DESCENT_WITHIN_30_BACKTRACKS'
        break
else:
    report['termination']='ITERATION_LIMIT'
r,j,n,fi,fo=p.evaluate(x)
report['final_norm']=n
report['final_margin_si']=float(p.margins(x)[margin_index])
report['final_mass_min_normalized']=x.reshape(-1,3)[:,:2].min(axis=0).tolist()
report['max_delta_over_old_scale']=float(np.max(np.abs(x-x0)/p.scale))
report['original_residual_pass']=bool(n<=1e-8)
report['calls']=p.calls
report['final_boundary_mass_flux_normalized']=[fi[:2].tolist(),fo[:2].tolist()]
suffix='' if args.fd_factor==1 else f'_fd{args.fd_factor:g}'
np.savez(ROOT/(args.case+'_active_solution'+suffix+'.npz'),x=x,r=r,scale=p.scale)
(ROOT/(args.case+'_active_report'+suffix+'.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('iterations','invalid_proposals')},indent=2))
for item in report['iterations']:
    print({k:v for k,v in item.items() if k!='gradient'})
