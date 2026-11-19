import argparse
import json
from offline_probe import Probe, InvalidState, ROOT, np

parser = argparse.ArgumentParser()
parser.add_argument('--case', default='D_SS1')
parser.add_argument('--mode', type=int, default=1)
args = parser.parse_args()
out = {'case': args.case, 'mode': args.mode}
try:
    p = Probe(args.case,args.mode)
    out['metadata'] = p.meta.tolist()
    r,j,n,_,_ = p.evaluate(p.trial)
    out['initial_norm'] = n
    saved = np.loadtxt(ROOT/args.case/'last_rejected_state.dat',skiprows=2)[:,-3:].ravel()
    out['saved_residual_max_difference'] = float(np.max(np.abs(r-saved)))
    out['margins'] = p.margins(p.trial).tolist()
    out['interface'] = p.interface(p.trial).tolist()
    dx = np.linalg.solve(j,-r)
    proposals=[]
    for k in (0,17):
        x = p.trial + 2.**-k*dx
        item={'lambda':2.**-k,'margins':p.margins(x).tolist(), 'mass_min':x.reshape(-1,3)[:,:2].min(axis=0).tolist()}
        try:
            item['norm']=p.evaluate(x,False)[2]
            item['status']=0
        except InvalidState as error:
            item.update(error.details)
        proposals.append(item)
    out['proposals']=proposals
    np.savez(ROOT/f'{args.case}_mode{args.mode}_initial.npz',x=p.trial,old=p.old,r=r,j=j,dx=dx,scale=p.scale)
except InvalidState as error:
    out['failure']=error.details
(ROOT/f'{args.case}_mode{args.mode}_check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
