#!/usr/bin/env python3
"""Consolidate existing final-state probe summaries; never run a transient."""
from __future__ import annotations
import json
from pathlib import Path
from build_final_domain_probe import read_metrics


def main() -> int:
    base = Path(__file__).resolve().parent
    root = base.parents[1]
    cases = []
    preceding = []
    probe_paths = list((base / 'final_ss3_cases').glob('*/summary.json'))
    probe_paths += list((base / 'final_osc_cases').glob('*/summary.json'))
    for path in sorted(probe_paths):
        summary = json.loads(path.read_text())
        metrics = summary.get('metrics', {})
        age = summary.get('snapshot_age_relative_to_run_end_s', 0)
        if abs(age) > 1e-10:
            preceding.append({'case': summary['case'], 'snapshot_time_s': metrics['time_s'],
                              'run_end_time_s': summary['run_status']['reached_time_s'],
                              'snapshot_age_s': age, 'classification': summary['classification'],
                              'summary': str(path.relative_to(root))})
            continue
        regime = summary['metadata']['regime'].removeprefix('D_')
        entry = {'case': summary['case'], 'regime': regime, 'eos': summary['metadata']['gas_eos'],
                 'pipe_cells': metrics['pipe_cells'], 'riser_cells': metrics['riser_cells'],
                 'dt_max_s': summary['metadata']['dt_max_s'],
                 'reached_time_s': summary['run_status']['reached_time_s'],
                 'requested_time_s': summary['run_status']['requested_time_s'],
                 'run_status': summary['run_status']['status'], 'classification': summary['classification'],
                 'diagnostic_input': summary['input_state_source'], 'point_failures': metrics['point_failures'],
                 'max_spectral_imag': metrics['max_spectral_imag'], 'min_ul_m_s': metrics['critical_ul_m_s'],
                 'min_ul_cell_extended': metrics['minimum_ul_cell_extended'],
                 'source_mismatches': summary['source_mismatches_against_run'],
                 'summary': str(path.relative_to(root))}
        if metrics.get('first_failed_face'):
            entry['first_failed_face'] = [metrics['first_failed_face'], metrics['first_failed_face'] + 1]
            entry['first_failed_face_status'] = metrics['first_failed_face_status']
        if summary['classification'] == 'MATCHED_JUNCTION_COCURRENT_TRACE_LIMIT':
            junction_path = base / 'final_junction_cases' / summary['case'] / 'summary.json'
            if junction_path.exists():
                junction = json.loads(junction_path.read_text())
                entry.update(junction_summary=str(junction_path.relative_to(root)),
                             junction_extended_cells=junction['junction_extended_cells'],
                             jg_trace_m_s=junction['profiles']['2']['jg_m_s'],
                             jg_newton_candidate_m_s=junction['profiles']['10']['jg_m_s'])
        cases.append(entry)
    for path in sorted((base / 'final_ss1_cases').glob('*/summary.json')):
        summary = json.loads(path.read_text())
        case = root / 'validation' / 'results' / summary['case']
        manifest = json.loads((case / 'manifest.json').read_text())
        status = read_metrics(case / 'run_status.txt')
        profile = summary['profiles']['2']
        candidate = summary['profiles'].get('10', {})
        classification = ('MATCHED_JUNCTION_COCURRENT_TRACE_LIMIT'
                          if candidate.get('jg_m_s', 0) < 0 and not summary['point_failures']
                          else 'JUNCTION_LIMIT_REQUIRES_DETAIL_REVIEW')
        if int(status['status']) == 3021:
            classification = 'MATCHED_JUNCTION_TRACE_SEARCH_LIMIT'
        cases.append({'case': summary['case'], 'regime': 'SS1', 'eos': manifest['controls']['eos'],
                      'pipe_cells': manifest['controls']['cells'], 'riser_cells': manifest['controls']['riser_cells'],
                      'dt_max_s': manifest['controls']['dt'], 'reached_time_s': status['reached_time_s'],
                      'requested_time_s': status['requested_time_s'], 'run_status': status['status'],
                      'classification': classification, 'diagnostic_input': 'rejected',
                      'snapshot_time_s': summary['time_s'],
                      'snapshot_age_relative_to_run_end_s': status['reached_time_s'] - summary['time_s'],
                      'point_failures': len(summary['point_failures']), 'minimum_cell_alpha': profile['minimum_cell_alpha'],
                      'junction_extended_cells': summary['junction_extended_cells'], 'jg_trace_m_s': profile['jg_m_s'],
                      'jg_newton_candidate_m_s': candidate.get('jg_m_s'), 'summary': str(path.relative_to(root))})
    cases.sort(key=lambda item: (item['regime'], item['eos'], item['pipe_cells'], -item['dt_max_s']))
    result = {'scope': 'Local diagnoses of final finite-volume runs, not experimental validation or convergence of slug cycles. Different run failures have distinct causes.',
              'cases': cases, 'older_rejected_snapshots_not_used_as_final_state': preceding,
              'trace_search_qualification': 'The steady alpha scan is limited to abs(logit(alpha))<=28. NO_ROOT status 21 does not prove that no positive root exists below its lower endpoint; the local Newton probes separately record proposed negative gas traces.',
              'trace_search_counterprobe': {'summary': 'validation/diagnostics/alpha_search_limit_counterprobe/summary.json',
                                           'scope': 'Synthetic local positive-gas trace, not an additional campaign case or the exact rejected Newton attempt.',
                                           'independent_alpha': 3.6117896990909259e-13,
                                           'production_alpha_min': 6.914400106935422e-13,
                                           'helper_status': 21,
                                           'closure_status_at_independent_root': 0},
              'source_provenance': 'Per-case build_provenance.json, snapshot_provenance.json and run manifest accompany the summaries.',
              'production_changes_from_diagnosis': False}
    (base / 'final_domain_summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    for item in cases:
        print(item['regime'], item['eos'], item['pipe_cells'], item['dt_max_s'], item['reached_time_s'],
              item['run_status'], item['classification'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
