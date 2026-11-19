from pathlib import Path
from collections import Counter
import argparse
import csv
import hashlib
import json

parser=argparse.ArgumentParser()
parser.add_argument('--label',default='active_set_smokes_v1')
args=parser.parse_args()
root=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119')/args.label
manifest=json.loads((root/'manifest.json').read_text())
suffixes=['l','g','l_pipe','g_pipe','l_riser','g_riser']
mass_columns=['Ml_kg','Mg_kg','Ml_pipe_kg','Mg_pipe_kg','Ml_riser_kg','Mg_riser_kg']
report={'label':args.label,'binary_sha256':manifest['binary_sha256'],'runs':[]}
for completed in manifest['runs']:
    case=completed['case']
    folder=root/case
    def complete_rows(name):
        rows=list(csv.DictReader((folder/name).open()))
        return [row for row in rows if None not in row and None not in row.values()]
    ledger=complete_rows('conservation.csv')
    steps=complete_rows('steps.csv')
    active=complete_rows('active_set.csv')
    proposals=complete_rows('rejected_proposals.csv')
    keys=['balance_'+name+'_kg' for name in suffixes]
    residual_keys=['mass_residual_sum_'+name+'_kg' for name in suffixes]
    final=ledger[-1]
    record={'case':case,'returncode':completed['returncode'],'timed_out':completed['timed_out'],
            'run_status':completed.get('run_status.txt'), 'last_rejected_iteration':completed.get('last_rejected_iteration.txt'),
            'snapshot_role':'last prior rejection, not a terminal state' if completed['timed_out'] else 'terminal linearization',
            'progress_role':'last completely flushed accepted records; lower bound at external timeout' if completed['timed_out'] else 'terminal run records',
            'balance_order':suffixes,'final_balance_kg':[float(final[k]) for k in keys],
            'max_abs_balance_kg':[max(abs(float(row[k])) for row in ledger) for k in keys],
            'max_balance_minus_integrated_residual_kg':[max(abs(float(row[k])-float(row[r])) for row in ledger) for k,r in zip(keys,residual_keys)],
            'min_inventory_kg':[min(float(row[k]) for row in ledger) for k in mass_columns],
            'max_junction_phase_flux_discrepancy_kg_s':[max(abs(float(row[f'mdot_{phase}_base_left_kg_s'])-float(row[f'mdot_{phase}_base_right_kg_s'])) for row in ledger) for phase in ('l','g')],
            'proposal_status_counts':dict(Counter(row['assembly_status'] for row in proposals)),
            'active_calls':len(active),'accepted_active_proposals':sum(row['accepted']=='T' for row in active),
            'first_active_record':active[0] if active else None,'last_active_record':active[-1] if active else None,
            'last_accepted_step':steps[-1] if steps else None,
            'last_proposals':proposals[-20:],
            'source_hashes':manifest['source_hashes'],
            'files_sha256':{name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ('last_rejected_state.dat','last_rejected_iteration.txt','run_status.txt','conservation.csv','steps.csv','active_set.csv','rejected_proposals.csv') if (folder/name).exists()}}
    report['runs'].append(record)
(root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
for row in report['runs']:
    print(row['case'],row['final_balance_kg'],row['max_balance_minus_integrated_residual_kg'],row['max_junction_phase_flux_discrepancy_kg_s'])
