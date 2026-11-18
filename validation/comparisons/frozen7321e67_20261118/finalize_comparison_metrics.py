#!/usr/bin/env python3
from pathlib import Path
import argparse,csv,subprocess,json,hashlib,sys
parser=argparse.ArgumentParser(description='Summarize the frozen comparison and retain failed-run classifications.')
parser.add_argument('--repo',type=Path,default=Path('/home/leonardoab/Documentos/PRIMS-F-comparison-20261118'))
root=parser.parse_args().repo.resolve()
output=root/'validation/comparisons/frozen7321e67_summary'
output.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'scripts'))
from analyze_validation import measure,read_metadata
from analyze_osc_envelopes import parse_status

folders=sorted((root/'validation/results').glob('comparison*_20261118_*'))
assert len(folders)==13,len(folders)
assert all(json.loads((p/'manifest.json').read_text())['status'] not in ['RUNNING'] for p in folders)
rows=[];inputs={}
physics_keys=['x','z','lp','beta_deg','d','eps','rho_l0','p_l0','cl','rg','mu_l','mu_g','mul','mug','sigma','g',
              'ps','jl','jg','tgas','gamma_g','gas_eos','gas_flow_basis','cg_fixed','p_standard','t_standard']
for folder in folders:
    metric,data=measure(folder,900.,100.)
    manifest=json.loads((folder/'manifest.json').read_text())
    status=parse_status((folder/'run_status.txt').read_text()) if (folder/'run_status.txt').exists() else {}
    case=manifest['controls']['case'];closure=manifest['controls']['closure']
    physics={k:metric['parameters'].get(k) for k in physics_keys}
    assert all(v is not None for v in physics.values()),physics
    if case in inputs:assert physics==inputs[case],(case,physics,inputs[case])
    else:inputs[case]=physics
    sha=manifest['source_sha256']
    assert all(hashlib.sha256((root/k).read_bytes()).hexdigest()==v for k,v in sha.items())
    row={'run':folder.name,'case':case,'closure':closure,'cells_each_section':manifest['controls']['cells'],
         'dt_max_s':manifest['controls']['dt'],'requested_time_s':manifest['controls']['time'],
         'status':manifest['status'],'solver_status':status.get('status'),
         'completed':metric['completed'],'reached_time_s':metric['reached_time_s'],
         'measurement_window_s':[metric['window_start_s'],metric['window_end_s']],
         'wall_seconds':manifest['wall_seconds'],'mass':metric['mass'],
         'sources_unchanged':manifest['sources_unchanged'],'binary_unchanged':manifest['binary_unchanged'],
         'monitor_sha256':hashlib.sha256((folder/'monitor.csv').read_bytes()).hexdigest(),
         'metrics_window_complete':metric['completed'] and metric['reached_time_s']>=2400}
    if row['metrics_window_complete']:
        row['cycles']=metric['cycles'];row['phase_indicators']=metric['phase_indicators']
    else:
        row['cycles']={'available':False,'reason':'2400 s and window 900-2400 not completed; no pilot-window substitution'}
    rows.append(row)
summary={'freeze_commit':'7321e6744932200176a0a6628b07e195d342071c','runs':rows,
         'same_physical_inputs_verified_within_case':True,'physical_inputs':inputs,
         'interpretation':'Numerical evidence; failed runs do not qualify a period or regime. Stable/OSC alone are not experimental validation.'}
original=next(p for p in folders if p.name.startswith('comparison_20261118_D_OSC_C3_'))
extension=next(p for p in folders if p.name.startswith('comparison_extended_20261118_D_OSC_C3_'))
old=list(csv.reader((original/'monitor.csv').open()));new=list(csv.reader((extension/'monitor.csv').open()))
limit=next(i for i,v in enumerate(old[1:],1) if float(v[0])>=2400)
assert old[:limit]==new[:limit]
summary['extension_reproduces_original_monitor']={'rows_including_header':limit,'bit_identical':True,
 'excluded_original_final_shortened_step_time_s':float(old[limit][0])}
(output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
for closure in [3,5]:
    osc=[p for p in folders if '_D_OSC_C'+str(closure)+'_' in p.name and not p.name.startswith('comparison_extended_')]
    cmd=['python3','scripts/analyze_osc_envelopes.py','--output',str(output/f'OSC_C{closure}')]
    for p in osc:cmd+=['--run-dir',str(p)]
    result=subprocess.run(cmd,cwd=root,text=True,capture_output=True)
    (output/f'OSC_C{closure}_analysis.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr
    cmd=['python3','scripts/analyze_phase_sequence.py','--output',str(output/f'phases_C{closure}'),
         '--window-start','900','--window-end','2400']
    for p in osc:cmd+=['--run-dir',str(p)]
    result=subprocess.run(cmd,cwd=root,text=True,capture_output=True)
    (output/f'phases_C{closure}_analysis.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr
extended=[p for p in folders if p.name.startswith('comparison_extended_')]
for p in extended:
    for script,out,args in [('analyze_osc_envelopes.py','OSC_C3_extended',[]),
                            ('analyze_phase_sequence.py','phases_C3_extended',['--window-start','3000','--window-end','3600'])]:
        cmd=['python3','scripts/'+script,'--run-dir',str(p),'--output',str(output/out),*args]
        result=subprocess.run(cmd,cwd=root,text=True,capture_output=True)
        (output/(out+'_analysis.log')).write_text(result.stdout+result.stderr)
        assert result.returncode==0,result.stderr
print(json.dumps({'runs':len(rows),'completed':sum(r['completed'] for r in rows),'output':str(output)},indent=2))
