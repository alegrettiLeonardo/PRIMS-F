#!/usr/bin/env python3
"""Pure-phase equilibrium at the outlet backflow limit, through the main CLI.

Zero feed, G=0, uniform state at the separator pressure and no external
reservoir composition (outlet_reservoir_phase=0). The equilibrium is admissible
and must evolve without step collapse; the effective backflow produced by an
unbalanced column (G>0, POINT source) must stay refused by status 934, never
converted to 935 or accepted with an invented composition.
"""
import argparse,csv,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'));sys.path.insert(0,str(ROOT/'scripts'))
from test_pipe_regime_cli import make_deck
from analyze_accepted_conservation import audit_ledger

def run(binary,out,name,extra):
 o=out/name;deck=out/(name+'.nml')
 params=dict(nodes=8,n_riser=8,jl=0,mul=0,jg=0,mug=0,G=0.0,time_dim=2.0,single_phase_mode="'REDUCED'",
             source_scheme="'POINT'",outlet_reservoir_phase=0,gas_eos="'ISOTHERMAL'",tol_res=1e-10,tol_inc=1e-10,
             newton_max=25,dt_start_s=.001,dt_max_s=.05,dt_min_s=1e-8,field_interval_s=.5,output_interval_s=.5,
             out_dir="'"+str(o.resolve())+"'")
 params.update(extra)
 deck.write_text(make_deck((ROOT/'input/benchmarks/D_Stable.nml').read_text(),params))
 r=subprocess.run([str(binary.resolve()),str(deck.resolve())],cwd=ROOT,capture_output=True,text=True,timeout=600,
                  env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
 (out/(name+'.log')).write_text(r.stdout+r.stderr)
 status=dict(l.split('=',1) for l in (o/'run_status.txt').read_text().split() if '=' in l)
 return r.returncode,status,o

def main():
 p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();a.output.mkdir(exist_ok=False);records=[]
 for phase,ini in [('liquid',"'LIQUID'"),('gas',"'GAS'")]:
  for scheme in ['POINT','HYDROSTATIC']:
   name=f'{phase}_{scheme}_equilibrium'
   rc,st,o=run(a.binary,a.output,name,dict(initial_phase=ini,source_scheme=f"'{scheme}'"))
   assert rc==0 and st['status']=='0' and st['termination_reason']=='COMPLETED',(name,st)
   assert int(st['rejected_steps'])==0,(name,st)
   ledger=audit_ledger(o/'conservation.csv');assert ledger['arithmetic_contract_passes'],name
   for key,v in ledger['balances'].items():
    assert v['max_abs_balance_kg']<=1e-12*max(1.0,v.get('initial_inventory_kg') or 1.0),(name,key,v)
   with(o/'boundary_conditions.csv').open()as f:rows=list(csv.DictReader(f))
   assert rows and all(x['status_in']==x['status_out']=='0' for x in rows),name
   records.append(dict(case=name,status=st['status'],accepted=st['accepted_steps'],rejected=st['rejected_steps']));print('PASS',name)
 # Effective physical backflow without external data stays refused.
 name='liquid_POINT_gravity_backflow'
 rc,st,o=run(a.binary,a.output,name,dict(initial_phase="'LIQUID'",G=9.81,time_dim=.05))
 assert rc!=0 and st['status']!='0',(name,st)
 with(o/'rejected_proposals.csv').open()as f:trials=list(csv.DictReader(f))
 codes={t['assembly_status'].strip() for t in trials}
 assert '934' in codes and '935' not in codes,(name,codes)
 records.append(dict(case=name,status=st['status'],termination=st['termination_reason'],proposal_codes=sorted(codes)));print('PASS',name)
 (a.output/'summary.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
