"""Condense the two frozen-baseline repeats and their three archived comparators.

Run after the standard envelope/phase/comparison analyzers. Full raw data stay
in their original directories; this summary is numerical evidence, not a regime
classification or experimental validation.
"""
import argparse
import hashlib
import importlib.metadata
import json
import platform
import statistics
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo', type=Path, required=True)
parser.add_argument('--analysis', type=Path, required=True)
args = parser.parse_args()
sys.path.insert(0, str(args.repo / 'scripts'))
from validation_evidence import load_monitor, time_statistics


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


env = json.loads((args.analysis / 'osc_envelopes.json').read_text())
for relative, expected in env['cohort']['source_sha256'].items():
    if digest(args.repo / relative) != expected:
        raise RuntimeError(f'Frozen source differs from run provenance: {relative}')
phase = {x['run']: x for x in json.loads((args.analysis / 'phase_sequence.json').read_text())['runs']}
comparison = json.loads((args.analysis / 'comparison.json').read_text())
metrics = {x['run']: x for x in comparison['runs']}
channels = {
    'p_junction_Pa': 'Pa', 'holdup_liquid_base': '1', 'holdup_liquid_out': '1',
    'mdot_l_base_kg_s': 'kg/s', 'mdot_g_base_kg_s': 'kg/s',
    'mdot_l_out_kg_s': 'kg/s', 'mdot_g_out_kg_s': 'kg/s',
    'Ml_kg': 'kg', 'Mg_kg': 'kg', 'Ml_riser_kg': 'kg', 'Mg_pipe_kg': 'kg',
    'Ml_pipe_kg': 'kg', 'Mg_riser_kg': 'kg',
}
cases = []
for case in env['cases']:
    directory = Path(case['source_directory'])
    manifest = json.loads((directory / 'manifest.json').read_text())
    if not case['completed'] or manifest['status'] != 'COMPLETED':
        raise RuntimeError(f'Refuse final summary of incomplete case: {directory}')
    data, warnings = load_monitor(directory / 'monitor.csv')
    final = case['last_ten_cycles']
    observation = {}
    if final.get('available'):
        bounds = [r['peak_time_s'] for r in final['peak_records']]
        for key, unit in channels.items():
            blocks = [time_statistics(data, key, a, b) for a, b in zip(bounds[:-1], bounds[1:])]
            if any(block is None or block['coverage_fraction'] < 1 - 1e-8 for block in blocks):
                observation[key] = {'available': False}
                continue
            means = [block['mean'] for block in blocks]
            excursions = [block['peak_to_peak'] for block in blocks]
            first_mean, last_mean = statistics.mean(means[:5]), statistics.mean(means[-5:])
            first_range, last_range = statistics.mean(excursions[:5]), statistics.mean(excursions[-5:])
            observation[key] = {
                'available': True, 'unit': unit,
                'first_five_mean': first_mean, 'last_five_mean': last_mean,
                'mean_relative_change': (last_mean - first_mean) / max(abs(first_mean), abs(last_mean), 1e-30),
                'first_five_mean_excursion': first_range, 'last_five_mean_excursion': last_range,
                'excursion_relative_change': (last_range - first_range) / max(abs(first_range), abs(last_range), 1e-30),
                'cycle_means': means, 'cycle_excursions': excursions,
            }
    p = phase[case['case']]
    item = {
        'name': case['case'], 'source_directory': case['source_directory'],
        'new_repeat': case['case'].startswith('repeat_20261118_frozen05ad799_'),
        'completed': True, 'observed_end_s': case['observed_end_s'],
        'nodes_each_section': case['nodes_pipe'], 'dt_max_s': case['dt_max_s'],
        'wall_seconds': manifest['wall_seconds'],
        'binary_sha256': manifest['binary_sha256'],
        'sources_unchanged': manifest['sources_unchanged'],
        'binary_unchanged': manifest['binary_unchanged'],
        'source_sha256': manifest['source_sha256'],
        'namelist_sha256': manifest['namelist_sha256'],
        'data_quality_warnings': warnings,
        'pressure_900_2400': case['pressure_900_2400'],
        'last_ten_pressure_cycles': {k: v for k, v in final.items() if k != 'peak_records'},
        'last_ten_phase_observables': observation,
        'mass_balance_full_run': metrics[case['case']]['mass'],
        'phase_window': p['window'], 'phase_classification': p['classification'],
        'phase_event_counts': {k: len(v) for k, v in p['events'].items()},
        'gas_base_observations': p.get('gas_base_observations'),
        'phase_signal_statistics': p['measured_signals'],
        'phase_thresholds': p.get('resolved_thresholds'),
        'raw_files': {name: {'sha256': digest(directory / name), 'bytes': (directory / name).stat().st_size}
                      for name in ('monitor.csv', 'fields.csv', 'manifest.json', 'run_metadata.txt', 'run_status.txt')
                      if (directory / name).is_file()},
    }
    cases.append(item)
summary = {
    'schema': 'PRIMS-F.frozen-baseline-repeat.v1',
    'summary_script_sha256': digest(Path(__file__)),
    'frozen_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.repo, text=True).strip(),
    'current_compiler': subprocess.check_output(['gfortran', '--version'], text=True).splitlines()[0],
    'current_python': platform.python_version(),
    'current_platform': platform.platform(),
    'current_analysis_packages': {name: importlib.metadata.version(name)
                                  for name in ('numpy', 'scipy', 'matplotlib')},
    'historical_environment': json.loads((args.repo / 'validation/environment.json').read_text()),
    'cohort_same_physics_and_source_hashes': True,
    'current_frozen_worktree_matches_every_recorded_source_hash': True,
    'cross_compiler_comparison': True,
    'restart_policy': 'Both repeats started from zero; archived interrupted runs were not concatenated.',
    'interpretation': 'Numerical refinement/phase observables only. No automatic physical regime or experimental validation.',
    'cases': cases,
    'comparisons': [{k: v for k, v in pair.items() if k not in ('coarse_summary', 'fine_summary')}
                    for pair in comparison['comparisons']],
    'analysis_file_sha256': {name: digest(args.analysis / name)
                            for name in ('osc_envelopes.json', 'phase_sequence.json', 'comparison.json')},
}
(args.analysis / 'frozen_repeat_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
print(json.dumps({'completed_cases': len(cases), 'new_repeats': sum(c['new_repeat'] for c in cases)}))
