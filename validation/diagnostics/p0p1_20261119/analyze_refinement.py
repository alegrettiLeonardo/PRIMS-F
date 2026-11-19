#!/usr/bin/env python3
"""End-cycle diagnostics for the isolated common-initial-value refinement."""
import argparse,csv,json,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
from analyze_osc_envelopes import last_ten_cycles,pressure_envelope
from analyze_accepted_conservation import audit_ledger
from cycle_waveform import wave_difference

def data(path):
 with path.open()as f:rows=list(csv.DictReader(f))
 # Keep the last representation if legacy finalization repeats a terminal time.
 rows=list({float(x['time_s']):x for x in rows}.values())
 return {k:np.array([float(x[k])for x in rows])for k in rows[0]}

def integral(t,y,l,r):
 inside=(t>l)&(t<r);x=np.r_[l,t[inside],r];v=np.r_[np.interp(l,t,y),y[inside],np.interp(r,t,y)]
 return float(np.trapz(v,x))

def cycle_markers(summary):
 edges=sorted(set([x['peak_time_s']for x in summary['peak_records']]+[summary['window_end_s']]))
 if len(edges)!=11 or any(r<=l for l,r in zip(edges[:-1],edges[1:])):
  raise ValueError('ten complete intervals require eleven distinct ordered markers')
 return edges

def summarize(folder):
 d=data(folder/'monitor.csv');t=d['time_s'];p=d['p_riser_base_Pa'];m=last_ten_cycles(t,p,pressure_envelope(t,p))
 result=dict(run=str(folder),time_s=float(t[-1]),cycles=m,settled=False)
 if not m['available']:return result
 edges=cycle_markers(m)
 periods=np.diff(edges);amp=np.array([x['semi_amplitude_Pa']for x in m['peak_records']if x['peak_time_s']<edges[-1]]);means=np.array([integral(t,p,l,r)/(r-l)for l,r in zip(edges[:-1],edges[1:])])
 variation={name:float(np.ptp(x)/abs(np.mean(x)))for name,x in [('period',periods),('amplitude',amp),('mean',means)]}
 result.update(settled=all(x<=.005 for x in variation.values()),settling_variation=variation,period_s=float(np.mean(periods)),semi_amplitude_Pa=float(np.mean(amp)),mean_Pa=float(np.mean(means)),window_s=[edges[0],edges[-1]])
 channels=['p_riser_base_Pa','p_junction_Pa','holdup_liquid_base','holdup_liquid_out','mdot_l_base_kg_s','mdot_g_base_kg_s','mdot_l_out_kg_s','mdot_g_out_kg_s','Ml_kg','Mg_kg','Ml_pipe_kg','Mg_pipe_kg','Ml_riser_kg','Mg_riser_kg']
 waves={}
 for key in channels:
  y=np.array([np.interp(l+(r-l)*np.arange(256)/256,t,d[key])for l,r in zip(edges[:-1],edges[1:])]);waves[key]=dict(mean=y.mean(0).tolist(),std_across_cycles=y.std(0).tolist())
 result['phase_waveforms']=waves
 ledger=audit_ledger(folder/'conservation.csv');result['balances']=ledger['balances'];result['arithmetic_contract_passes']=ledger['arithmetic_contract_passes'];result['junction_mismatch_kg']=ledger['max_junction_integral_mismatch_kg']
 c=data(folder/'conservation.csv');result['per_cycle_balances_kg']={}
 for phase in ['l','g']:
  for region,mass,fi,fo in [('system','',f'M{phase}_in_integral_kg',f'M{phase}_out_integral_kg'),('pipe','_pipe',f'M{phase}_in_integral_kg',f'M{phase}_base_left_integral_kg'),('riser','_riser',f'M{phase}_base_right_integral_kg',f'M{phase}_out_integral_kg')]:
   masskey=f'M{phase}{mass}_kg';balance=c[masskey]-c[fi]+c[fo]
   result['per_cycle_balances_kg'][phase+'_'+region]=np.diff(np.interp(edges,c['time_s'],balance)).tolist()
 result['cycle_balance_definition']='linear interpolation of accepted BE inventories and cumulative accepted flux integrals at the same peak markers; no re-integration of sampled monitors'
 steps=data(folder/'steps.csv');result['actual_dt_s']=dict(min=float(np.min(steps['dt_s'][1:])),max=float(np.max(steps['dt_s'])),mean=float(np.mean(steps['dt_s'][1:])),accepted=len(steps['dt_s'])-1)
 result['hashes']={name:hashlib.sha256((folder/name).read_bytes()).hexdigest()for name in ['monitor.csv','steps.csv','conservation.csv','checkpoint.dat']}
 return result

def main():
 p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,required=True);p.add_argument('--extensions',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
 result=[]
 for folder in sorted(a.runs.glob('N*')):
  if not(folder/'run_status.txt').exists():continue
  use=folder
  if a.extensions and(a.extensions/folder.name/'run_status.txt').exists():use=a.extensions/folder.name
  d=summarize(use);d['base_run']=str(folder);n,dt=folder.name.split('_dt');d.update(n=int(n[1:]),dt_max_s=float(dt));result.append(d)
 comparisons=[]
 for axis,cohort in [('space',sorted([x for x in result if x['dt_max_s']==.125],key=lambda x:x['n'])),('time',sorted([x for x in result if x['n']==51],key=lambda x:-x['dt_max_s']))]:
  for coarse,fine in zip(cohort,cohort[1:]):
   if 'period_s'not in coarse or'period_s'not in fine:continue
   d=dict(axis=axis,coarse=Path(coarse['base_run']).name,fine=Path(fine['base_run']).name,both_settled=coarse['settled']and fine['settled'])
   for key in ['period_s','semi_amplitude_Pa','mean_Pa']:
    d[key+'_relative_difference']=abs(fine[key]-coarse[key])/max(abs(fine[key]),abs(coarse[key]))
   d['waveforms']={k:wave_difference(coarse['phase_waveforms'][k],fine['phase_waveforms'][k],fine['period_s'])for k in fine['phase_waveforms']}
   comparisons.append(d)
 fig,axes=plt.subplots(3,2,figsize=(10,9))
 for column,(axis,key,label) in enumerate([('space','n','space'),('time','dt_max_s','time')]):
  cohort=sorted([x for x in result if(x['dt_max_s']==.125 if axis=='space' else x['n']==51)],key=lambda x:x[key])
  for j,observable in enumerate(['period_s','semi_amplitude_Pa','mean_Pa']):
   good=[x for x in cohort if observable in x];ax=axes[j,column];xx=[57.4/x['n'] if axis=='space'else x[key]for x in good];ax.plot(xx,[x[observable]for x in good],'o-');ax.set_xlabel('pipe cell width (m)' if axis=='space'else 'dt_max (s)');ax.set_title(observable);ax.set_xscale('log',base=2);ax.invert_xaxis();ax.grid(alpha=.2)
 fig.tight_layout();fig.savefig(a.output/'refinement.png',dpi=160);plt.close(fig)
 orders=[]
 for axis in ['space','time']:
  cohort=sorted([x for x in result if(x['dt_max_s']==.125 if axis=='space'else x['n']==51)],key=lambda x:-57.4/x['n']if axis=='space'else-x['dt_max_s'])[-3:]
  if len(cohort)<3 or any(not x['settled']for x in cohort):continue
  h=np.array([57.4/x['n']if axis=='space'else x['dt_max_s']for x in cohort])
  for key in ['period_s','semi_amplitude_Pa']:
   y=np.array([x[key]for x in cohort]);ratio=(y[2]-y[1])/(y[1]-y[0]);entry=dict(axis=axis,observable=key,h=h.tolist(),values=y.tolist(),asymptotic_region_proven=False)
   try:
    def fun(p):return(h[2]**p-h[1]**p)/(h[1]**p-h[0]**p)-ratio
    order=brentq(fun,.05,8);r=h[1]/h[2];ext=y[2]+(y[2]-y[1])/(r**order-1);entry.update(observed_order=order,formal_extrapolate=float(ext),estimated_finest_error_relative=float(abs(ext-y[2])/abs(ext)))
   except(ValueError,ZeroDivisionError):entry['reason']='No positive observed order in [0.05,8] for the signed differences'
   orders.append(entry)
 (a.output/'summary.json').write_text(json.dumps(dict(runs=result,comparisons=comparisons,orders=orders,interpretation='new common-initial-value campaign; settling and discretization convergence are separate; no experimental validation or GCI claimed'),indent=2,allow_nan=False)+'\n')
 for x in result:print(x['n'],x['dt_max_s'],x.get('period_s'),x.get('semi_amplitude_Pa'),x.get('settling_variation'),x['settled'])
if __name__=='__main__':main()
