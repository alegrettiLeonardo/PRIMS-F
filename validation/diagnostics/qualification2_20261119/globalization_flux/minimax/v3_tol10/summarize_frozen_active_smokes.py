"""Compact results for independent frozen smoke manifests; never infer terminal status."""
from pathlib import Path
import argparse
import csv
import hashlib
import json


def rows(path):
    if not path.exists():
        return []
    with path.open() as stream:
        complete = (line for line in stream if line.endswith('\n'))
        return [row for row in csv.DictReader(complete)
                if None not in row and None not in row.values()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    report = {'schema': 'primsf_frozen_active_smoke_summary_v1', 'runs': []}
    suffixes = ('l', 'g', 'l_pipe', 'g_pipe', 'l_riser', 'g_riser')
    inventory = ('Ml_kg', 'Mg_kg', 'Ml_pipe_kg', 'Mg_pipe_kg', 'Ml_riser_kg', 'Mg_riser_kg')
    for case in ('D_SS1', 'D_SS3'):
        manifest_file = args.root / case / 'manifest.json'
        if not manifest_file.exists():
            continue
        manifest = json.loads(manifest_file.read_text())
        output = Path(manifest['output_directory'])
        steps = rows(output / 'steps.csv')
        conservation = rows(output / 'conservation.csv')
        active = rows(output / 'active_set.csv')
        rejected = rows(output / 'rejected_proposals.csv')
        item = {name: manifest.get(name) for name in (
            'case', 'source_commit', 'binary_sha256', 'state', 'termination_reason',
            'returncode', 'timed_out', 'elapsed_wall_s', 'child_cpu_seconds',
            'run_status', 'snapshot_role', 'progress_role')}
        item['manifest_sha256'] = hashlib.sha256(manifest_file.read_bytes()).hexdigest()
        item['last_complete_step'] = steps[-1] if steps else None
        item['balance_order'] = suffixes
        if conservation:
            item['final_conservation_position'] = {
                key: conservation[-1][key] for key in ('step', 'time_s', 'dt_s', 'residual_norm')}
            item['buffer_note'] = 'Each stream is truncated to its last newline-terminated complete row; steps and conservation may end at different accepted steps after external timeout.'
            item['final_balance_kg'] = [float(conservation[-1]['balance_' + phase + '_kg'])
                                        for phase in suffixes]
            item['max_balance_minus_integrated_residual_kg'] = [
                max(abs(float(row['balance_' + phase + '_kg']) -
                        float(row['mass_residual_sum_' + phase + '_kg']))
                    for row in conservation) for phase in suffixes]
            item['min_inventory_kg'] = [min(float(row[key]) for row in conservation)
                                        for key in inventory]
            item['max_junction_phase_flux_discrepancy_kg_s'] = [
                max(abs(float(row['mdot_' + phase + '_base_left_kg_s']) -
                        float(row['mdot_' + phase + '_base_right_kg_s']))
                    for row in conservation) for phase in ('l', 'g')]
        item['active_records'] = len(active)
        item['accepted_active_proposals'] = sum(row['accepted'] == 'T' for row in active)
        item['first_active_record'] = active[0] if active else None
        item['last_active_record'] = active[-1] if active else None
        item['last_rejected_proposals'] = rejected[-20:]
        if steps:
            item['last_50_accepted_steps'] = steps[-50:]
            item['max_accepted_residual_norm'] = max(float(row['residual_norm']) for row in steps)
        report['runs'].append(item)
    target = args.root / 'summary.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'summary': str(target), 'cases': len(report['runs'])}))


if __name__ == '__main__':
    main()
