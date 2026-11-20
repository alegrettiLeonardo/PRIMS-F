#!/usr/bin/env python3
"""Derive an auditable chronological analysis dataset from exact checkpoints.

No solver run_status is synthesized. Raw attempts, including discarded tails,
remain untouched. Integrals are independently recomputed from accepted BE
fluxes, never corrected to match inventories.
"""
import argparse,csv,json,sys,hashlib,shutil
import numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from analyze_accepted_conservation import Accumulator,audit_ledger
sys.path.insert(0,str(Path(__file__).parent))
from step_metadata import accepted_rows
CHANNELS=['p_riser_base_Pa','p_junction_Pa','holdup_liquid_base','holdup_liquid_out','mdot_l_base_kg_s','mdot_g_base_kg_s','mdot_l_out_kg_s','mdot_g_out_kg_s','Ml_kg','Mg_kg','Ml_pipe_kg','Mg_pipe_kg','Ml_riser_kg','Mg_riser_kg']

def rows(p):
 with p.open()as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def validate_operators(segments):
 operators=[]
 for segment in segments:
  manifest=json.loads((Path(segment['run'])/'manifest.json').read_text())
  operator=manifest.get('operator')
  if not isinstance(operator,str) or len(operator)<7:raise ValueError('A source operator SHA is required')
  operators.append(operator[:7])
 if len(set(operators))!=1:raise ValueError('Different operator versions cannot be spliced')
 return operators[0]

def compose(segments,out):
 validate_operators(segments)
 out.mkdir(exist_ok=False);combined=[];step_rows=[];mon=[];field_sources=[];provenance=[];count=0;counter_offsets={};previous=None
 for index,segment in enumerate(segments):
  folder=Path(segment['run']);offset=segment.get('offset_s',0.);last_step=segment.get('through_step');ledger=rows(folder/'conservation.csv')
  if last_step is not None:ledger=[x for x in ledger if int(x['step'])<=last_step]
  assert ledger and int(ledger[-1]['step'])==(last_step if last_step is not None else int(ledger[-1]['step']))
  if previous is not None:
   assert abs(float(ledger[0]['time_s'])+offset-float(previous['time_s']))<1e-10
   for phase in ['l','g']:
    for region in ['','_pipe','_riser']:
     k=f'M{phase}{region}_kg';error=abs(float(ledger[0][k])-float(previous[k]));assert error<=1e-11*max(abs(float(previous[k])),1e-20),('checkpoint inventory discontinuity',k,error)
  if index:
   ck=(Path(segments[index-1]['run'])/'checkpoint.dat').read_text().splitlines();hd=ck[1].split();nd=int(hd[1]);n=int(hd[0]);expected=np.array([list(map(float,x.split()[4:7]))for x in ck[2:]])[nd:n-nd]
   actual=np.loadtxt(folder/'initial_U.txt');assert np.array_equal(actual,expected), 'restart must use the exact exported conservative state, not only the same global inventory'
  local_end=float(ledger[-1]['time_s']);local_steps=int(ledger[-1]['step']);source_steps=rows(folder/'steps.csv');source_steps=[x for x in source_steps if int(x['step'])<=local_steps]
  for row in ledger[(1 if index else 0):]:
   row=dict(row);row['time_s']=format(float(row['time_s'])+offset,'.17e');row['step']=str(count+int(row['step']));combined.append(row)
  for row in accepted_rows(source_steps):
   row=dict(row);row['time_s']=format(float(row['time_s'])+offset,'.17e');row['step']=str(count+int(row['step']))
   for k in row:
    if k.endswith('_total'):row[k]=str(int(float(row[k])+counter_offsets.get(k,0)))
    elif k in ['cpu_total_s','wall_total_s']:row[k]=format(float(row[k])+counter_offsets.get(k,0),'.17e')
   step_rows.append(row)
  if source_steps:
   for k in source_steps[-1]:
    if k.endswith('_total')or k in ['cpu_total_s','wall_total_s']:counter_offsets[k]=counter_offsets.get(k,0)+float(source_steps[-1][k])
  for row in rows(folder/'monitor.csv'):
   if float(row['time_s'])>local_end:continue
   mon.append({k:(format(float(row[k])+offset,'.17e')if k=='time_s'else row[k])for k in ['time_s','dt_s',*CHANNELS]})
  field_sources.append((folder,offset,local_end));provenance.append(dict(run=str(folder),operator=json.loads((folder/'manifest.json').read_text()).get('operator'),source_manifest_sha256=sha(folder/'manifest.json'),offset_s=offset,through_step=local_steps,through_time_s=local_end,source_sha256={name:sha(folder/name)for name in ['conservation.csv','monitor.csv','steps.csv','fields.csv','checkpoint.dat']}))
  previous=combined[-1];count+=local_steps
 # Same accepted weights as the updates; different summation order is explicit.
 signed={(p,f):Accumulator()for p in ['l','g']for f in ['in','out','base_left','base_right']};absolute={k:Accumulator()for k in signed};residual={(p,r):Accumulator()for p in ['l','g']for r in ['','_pipe','_riser']};initial=combined[0]
 for i,row in enumerate(combined):
  dt=float(row['dt_s'])
  for phase in ['l','g']:
   for face in ['in','out','base_left','base_right']:
    if i:
     flux=float(row[f'mdot_{phase}_{face}_kg_s']);signed[phase,face].add(dt*flux);absolute[phase,face].add(dt*abs(flux))
    row[f'M{phase}_{face}_integral_kg']=format(signed[phase,face].value,'.17e');row[f'M{phase}_{face}_abs_integral_kg']=format(absolute[phase,face].value,'.17e')
   for region,inlet,outlet in [('', 'in','out'),('_pipe','in','base_left'),('_riser','base_right','out')]:
    if i:residual[phase,region].add(float(row[f'mass_residual_{phase}{region}_kg']))
    row[f'mass_residual_sum_{phase}{region}_kg']=format(residual[phase,region].value,'.17e')
    k=f'M{phase}{region}_kg';balance=float(row[k])-float(initial[k])-signed[phase,inlet].value+signed[phase,outlet].value
    row[f'balance_{phase}{region}_kg']=format(balance,'.17e')
 def write(name,records):
  with(out/name).open('w')as f:
   w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
 write('conservation.csv',combined);write('steps.csv',step_rows)
 mon=list({float(x['time_s']):x for x in mon}.values());write('monitor.csv',mon)
 with(out/'fields.csv').open('w')as output:
  writer=None;last_time=None
  for folder,offset,local_end in field_sources:
   with(folder/'fields.csv').open()as source:
    reader=csv.DictReader(source)
    if writer is None:writer=csv.DictWriter(output,fieldnames=reader.fieldnames);writer.writeheader()
    for row in reader:
     if float(row['time_s'])>local_end:continue
     row['time_s']=format(float(row['time_s'])+offset,'.17e');writer.writerow(row)
 shutil.copy2(Path(segments[-1]['run'])/'checkpoint.dat',out/'checkpoint.dat')
 audit=audit_ledger(out/'conservation.csv');assert audit['arithmetic_contract_passes'],audit['mismatches'][:3]
 manifest=dict(kind='derived analysis chronology; not a solver output/status',completed_analysis_history=True,segments=provenance,accepted_steps=count,reached_time_s=float(combined[-1]['time_s']),sum_rule='compensated independent re-integration of the original accepted BE fluxes and mass residuals; states are unchanged; step_balance remains the original measured quantity',balances=audit['balances'],arithmetic_contract_passes=audit['arithmetic_contract_passes'],files={p.name:sha(p)for p in out.iterdir()if p.is_file()})
 (out/'analysis_history.json').write_text(json.dumps(manifest,indent=2)+'\n');return manifest

def main():
 p=argparse.ArgumentParser();p.add_argument('--segments',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();m=compose(json.loads(a.segments.read_text()),a.output);print(m['accepted_steps'],m['reached_time_s'],m['arithmetic_contract_passes'])
if __name__=='__main__':main()
