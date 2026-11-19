from pathlib import Path
import argparse,csv,hashlib,json,math,sys
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_local import scalar
from audit_paths import ROOT
from audit_paths import OUT
sys.path.insert(0,str(ROOT/'scripts'))
from analyze_accepted_conservation import audit_ledger,audit_steps
def csvrows(path):
 with path.open() as f:return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
def normalized(error,denominator):
 return None if denominator==0 else error/denominator
def analyze_run(path):
 rows=csvrows(path/'boundary_contract.csv');accepted=rows[1:];mass=csvrows(path/'conservation.csv')
 assert len(accepted)==len(mass)-1
 for b,m in zip(rows,mass):assert b['step']==m['step'] and b['time_s']==m['time_s']
 data={'phase':{},'pressure_max_abs_Pa':max(abs(x['rp_Pa']) for x in rows),
  'flux_recompute_max_nd':max(x['flux_recompute_max_nd'] for x in rows),'accepted_steps':len(accepted),
  'final_time_s':rows[-1]['time_s'],'dt_min_used_s':min(x['dt_s'] for x in accepted),'dt_max_used_s':max(x['dt_s'] for x in accepted)}
 for phase in ['l','g']:
  errors=[x[f'rm_{phase}_kg_s'] for x in rows]
  targets=[x['mdot_l_Qtarget_kg_s' if phase=='l' else 'mdot_g_target_kg_s'] for x in rows]
  rel=[normalized(abs(e),abs(t)) for e,t in zip(errors,targets)]
  integral=math.fsum(x[f'rm_{phase}_kg_s']*x['dt_s'] for x in accepted)
  targetint=math.fsum(x['mdot_l_Qtarget_kg_s' if phase=='l' else 'mdot_g_target_kg_s']*x['dt_s'] for x in accepted)
  data['phase'][phase]={'max_instantaneous_abs_kg_s':max(map(abs,errors)),
   'max_instantaneous_relative':max((r for r in rel if r is not None),default=None),
   'integrated_signed_error_kg':integral,'relative_to_prescribed_mass':normalized(integral,abs(targetint)),
   'relative_to_initial_inventory':normalized(integral,mass[0][f'M{phase}_kg']),
   'prescribed_mass_kg':targetint,'initial_inventory_kg':mass[0][f'M{phase}_kg'],
   'inventory_change_kg':mass[-1][f'M{phase}_kg']-mass[0][f'M{phase}_kg']}
 ledger=audit_ledger(path/'conservation.csv');steps=audit_steps(path/'steps.csv',ledger)
 data['ledger']=ledger;data['steps']=steps
 for p in ['l','g']:
  data['phase'][p]['balances']={}
  for region in ['','_pipe','_riser']:
   key='global' if not region else region[1:];m0=mass[0][f'M{p}{region}_kg']
   proc=mass[-1][f'M{p}_in_abs_integral_kg']+mass[-1][f'M{p}_out_abs_integral_kg'] if not region else (mass[-1][f'M{p}_in_abs_integral_kg']+mass[-1][f'M{p}_base_left_abs_integral_kg'] if region=='_pipe' else mass[-1][f'M{p}_base_right_abs_integral_kg']+mass[-1][f'M{p}_out_abs_integral_kg'])
   e=mass[-1][f'balance_{p}{region}_kg']
   data['phase'][p]['balances'][key]={'error_kg':e,'relative_M0':normalized(abs(e),m0),'relative_processed':normalized(abs(e),proc),
     'max_step_defect_kg':max(abs(x[f'step_balance_{p}{region}_kg']) for x in mass),
     'junction_max_abs_mismatch_kg_s':max(abs(x[f'mdot_{p}_base_left_kg_s']-x[f'mdot_{p}_base_right_kg_s']) for x in mass)}
 return data
def main():
 result={'normalization_at_zero':{'zero_error_zero_target':normalized(0,0),'nonzero_error_zero_target':normalized(1,0),'absolute_errors_reported':True},'runs':{}}
 for path in sorted((OUT/'controls').glob('*'))+sorted((OUT/'windows').glob('*')):
  if path.is_dir() and (path/'boundary_contract.csv').exists():result['runs'][str(path.relative_to(OUT))]=analyze_run(path)
 for variant in ['standard','volume']:
  path=OUT/'windows'/variant;rows=csvrows(path/'monitor.csv');p=np.array([r['p_riser_base_Pa'] for r in rows]);t=np.array([r['time_s'] for r in rows])
  peaks,_=find_peaks(p,prominence=500,distance=100)
  result['runs']['windows/'+variant]['pressure_window']={'definition':'diagnostic 120-s reconstructed window, no regime/convergence claim','mean_time_trapezoid_Pa':np.trapezoid(p,t)/t[-1] if hasattr(np,'trapezoid') else np.trapz(p,t)/t[-1],
    'semi_excursion_Pa':.5*(max(p)-min(p)),'peak_times_s':t[peaks],'mean_peak_spacing_s':float(np.mean(np.diff(t[peaks]))),'complete_peak_intervals':len(peaks)-1}
 fig,ax=plt.subplots(3,1,figsize=(10,8),sharex=True)
 for v in ['standard','volume']:
  rows=csvrows(OUT/'windows'/v/'monitor.csv');t=[r['time_s'] for r in rows]
  ax[0].plot(t,np.array([r['p_riser_base_Pa'] for r in rows])/1e5,label=v)
  ax[1].plot(t,[r['Ml_kg'] for r in rows],label=v)
  rows=csvrows(OUT/'windows'/v/'boundary_contract.csv');t=[r['time_s'] for r in rows]
  ax[2].plot(t,[1e6*r['rm_l_kg_s'] for r in rows],label=v)
 ax[0].set_ylabel('p base (bar abs)');ax[1].set_ylabel('M líquido (kg)');ax[2].set_ylabel('r líquido (mg/s)');ax[2].set_xlabel('t da janela (s)')
 for axis in ax:axis.legend();axis.grid(alpha=.25)
 fig.suptitle('Janela OSC: mesmo estado inicial reconstruído; apenas contrato líquido alterado')
 fig.tight_layout();fig.savefig(OUT/'osc_boundary_comparison.png',dpi=160);plt.close(fig)
 # The single-cell hydrostatic reference is independent of the executed lower-face pressure.
 hydro=json.loads((OUT/'hydrostatic.json').read_text());a=__import__('audit_local').Audit(OUT/'release/audit.so');checks=[]
 for r in hydro:
  pc=r['executed_lower_face_pressure_Pa'];pb=r['hydrostatic_lower_face_reference_Pa'];theta=.7;dt=.001;rho=997.05+(pc-101325)/1498**2
  smid=-rho/a.par[1]*9.81*a.par[4]*math.sin(theta)
  ref=dt*((a.par[16]*a.par[2]-pb)/(a.par[1]*r['ds_nd'])-smid)
  executed=dt*((a.par[16]*a.par[2]-pc)/(a.par[1]*r['ds_nd'])-smid)
  checks.append({'ds_nd':r['ds_nd'],'independent_midpoint_hydrostatic_residual_nd':ref,'independent_executed_residual_nd':executed,'difference_from_Fortran':executed-r['momentum_residual_nd']})
 result['hydrostatic_independent']=checks
 (OUT/'analysis.json').write_text(json.dumps(scalar(result),indent=2,allow_nan=False)+'\n')
 for name,r in result['runs'].items():
  print(name,'max_mass_rel',max(x['max_instantaneous_relative'] or 0 for x in r['phase'].values()),'max_balance_M0',max(x['relative_M0'] or 0 for p in r['phase'].values() for x in p['balances'].values()),r.get('pressure_window',''))
if __name__=='__main__':main()
