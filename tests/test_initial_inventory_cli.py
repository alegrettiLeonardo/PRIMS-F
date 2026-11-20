#!/usr/bin/env python3
"""Explicit pure initial inventory versus zero-flow ambiguity, through main CLI."""
import argparse,csv,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'));sys.path.insert(0,str(ROOT/'scripts'))
from test_pipe_regime_cli import make_deck
from analyze_accepted_conservation import audit_ledger
def main():
 p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
 cases=[('gas_auto',{},False),('gas_declared',dict(initial_phase="'GAS'"),True),('liquid_declared',dict(initial_phase="'LIQUID'",jl=.001,mul=.99705,jg=0,mug=0,outlet_reservoir_phase=1),True),('incompatible',dict(initial_phase="'GAS'",mul=1),False),('strict_rejected',dict(initial_phase="'GAS'",single_phase_mode="'STRICT'"),False)]
 records=[]
 for name,extra,success in cases:
  out=a.output/name;deck=a.output/(name+'.nml');params=dict(nodes=8,n_riser=8,jl=0,mul=0,jg=.1,mug=.120403094823,time_dim=.1,single_phase_mode="'REDUCED'",source_scheme="'HYDROSTATIC'",outlet_reservoir_phase=2,gas_eos="'ISOTHERMAL'",tol_res=1e-10,tol_inc=1e-10,newton_max=25,dt_start_s=.001,dt_max_s=.01,dt_min_s=1e-8,field_interval_s=.01,output_interval_s=.01,out_dir="'"+str(out.resolve())+"'");params.update(extra)
  deck.write_text(make_deck((ROOT/'input/benchmarks/D_Stable.nml').read_text(),params));run=subprocess.run([str(a.binary.resolve()),str(deck.resolve())],cwd=ROOT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True,timeout=60);(a.output/(name+'.log')).write_text(run.stdout+run.stderr)
  assert (run.returncode==0)==success,name
  if success:
   ledger=audit_ledger(out/'conservation.csv');assert ledger['arithmetic_contract_passes']
   assert all((v['relative_inventory']is None and v['max_abs_balance_kg']==0)or(v['relative_inventory']is not None and v['relative_inventory']<=1e-6)for v in ledger['balances'].values())
   with(out/'boundary_conditions.csv').open()as f:rows=list(csv.DictReader(f))
   phase='1'if name.startswith('liquid')else'2';assert all(x['reduced_phase']==phase and x['nwaves_in']==x['nwaves_out']=='2'for x in rows)
   assert all(x['status_in']==x['status_out']=='0'for x in rows)
  records.append(dict(case=name,returncode=run.returncode,passed=True));print('PASS',name)
 (a.output/'summary.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
