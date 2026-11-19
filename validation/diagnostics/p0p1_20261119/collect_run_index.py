#!/usr/bin/env python3
"""Read-only delivery index; external timeout is never a physical failure."""
import argparse,csv,json,hashlib
from step_metadata import statistics
from pathlib import Path
GROUPS=['feed_qualified_1365094','p0_campaign_18aa7f7','guard_qualification_994736f','control_eos_tol_994736f','local_volume_controls_994736f','osc_feed_ab_1603419','pure_cli_71a0693','restart_refinement_18aa7f7','refinement_extensions_18aa7f7_r2','refinement_resumes_18aa7f7','controller_preserving_18aa7f7']
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();records=[]
 for group in GROUPS:
  for folder in sorted((a.results/group).glob('*')):
   if not folder.is_dir()or not(folder/'steps.csv').exists():continue
   manifest=json.loads((folder/'manifest.json').read_text())if(folder/'manifest.json').exists()else{}
   operator=manifest.get('operator',manifest.get('operator_sha'));operator_sources=[]
   if not operator:
    for name in ['manifest.json','protocol.json','PROTOCOL.json','summary.json']:
     source=folder.parent/name
     if not source.exists():continue
     parent=json.loads(source.read_text())
     if isinstance(parent,list):
      matching=[x for x in parent if isinstance(x,dict)and x.get('phase')==folder.name]
      candidate=matching[0].get('operator')if len(matching)==1 else None
     else:
      candidate=parent.get('operator')or parent.get('operator_sha')
      if not candidate and isinstance(parent.get('protocol'),dict):candidate=parent['protocol'].get('operator')
     if candidate:operator=candidate;operator_sources.append(dict(file=str(source),sha256=sha(source)));break
   with(folder/'steps.csv').open()as f:rows=list(csv.DictReader(f))
   if not rows:continue
   last=rows[-1];dt_stats=statistics(rows)
   status={k.strip():v.strip()for k,v in(x.split('=',1)for x in(folder/'run_status.txt').read_text().splitlines()if'='in x)}if(folder/'run_status.txt').exists()else None
   record=dict(run=str(folder),group=group,operator=operator,operator_parent_sources=operator_sources,external_timeout=manifest.get('timeout',False),solver_status=status,confirmed_last_accepted_time_s=float(last['time_s']),confirmed_last_accepted_step=int(last['step']),last_recorded_rejections=int(last['rejected_steps_total']),confirmed_cpu_lower_bound_s=float(last['cpu_total_s']),external_wall_seconds=manifest.get('wall_s',manifest.get('wall_seconds')),actual_dt_s=dt_stats,sha256={name:sha(folder/name)for name in ['manifest.json','input.nml','initial_U.txt','run_metadata.txt','run_status.txt','steps.csv','conservation.csv','boundary_conditions.csv','checkpoint.dat','monitor.csv','fields.csv']if(folder/name).exists()})
   if manifest.get('timeout'):
    record['interpretation']='Last recorded accepted data are a confirmed lower bound; no final solver status or physical failure is inferred. Raw files can have different flushed lengths.'
   records.append(record)
 build_receipts=[]
 for source in sorted(a.results.glob('frozen_*/build_manifest.json')):
  build=json.loads(source.read_text());build_receipts.append(dict(file=str(source),sha256=sha(source),operator=build.get('operator',build.get('sha')),binaries=build.get('binaries'),sources=build.get('sources'),compiler=build.get('compiler'),flags_reference='Makefile text in original build receipt; linked-library hashes in ENVIRONMENT.json'))
 result=dict(schema='PRIMS-F.P0P1.delivery.v1',build_receipts=build_receipts,records=records,limits=['Each original source/binary/input manifest retains its own version.','Analysis composites are indexed separately and are never labeled as original solver outputs.','CPU lower bounds for truncated attempts exclude any unflushed work; all paid attempts remain listed.'])
 with a.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('indexed',len(records),'runs')
if __name__=='__main__':main()
