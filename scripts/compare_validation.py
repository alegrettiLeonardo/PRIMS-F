#!/usr/bin/env python3
"""Compare independent mesh/time changes without declaring convergence from mass.

Examples:
  python scripts/compare_validation.py --metrics validation/analysis_final/metrics.json
  python scripts/compare_validation.py --results validation/results --pattern 'final_*'
  python scripts/compare_validation.py --pattern 'pilot5_*' --window-start 0 --output validation/diagnostics/pilot_comparison

The time control is dt_max_s, an adaptive-step ceiling, not a fixed time step.
Comparisons recompute monitor statistics on a common physical-time window.
Absent raw signals and incomplete runs remain explicit limitations.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from validation_evidence import (ROOT, analysis_window, finite, fmt, json_ready, load_runs,
                                 pressure_cycles, pressure_key, time_statistics, write_json)

# Exclude mesh counts and dt_max_s so each can be varied independently. The
# remaining physical/model/numerical choices must agree within a study family.
MATCH_KEYS = ('x', 'z', 'lp', 'beta_deg', 'd', 'eps', 'rho_l0', 'p_l0', 'cl', 'cg_fixed',
              'rg', 'mu_l', 'mu_g', 'sigma', 'g', 'ps', 'jl', 'jg', 'tgas', 'mul', 'mug',
              'gas_eos', 'gas_flow_basis', 'p_standard', 't_standard', 'closure_riser',
              'junction_model', 'dissipation_mode', 'flux_scheme', 'perturbation_amplitude',
              'mesh_ratio', 'use_muscl', 'quad_points', 'tol', 'tola', 'tol_res', 'tol_inc',
              'dt_start_s', 'dt_min_s', 'dt_growth', 'dt_shrink', 'dt_fac_min', 'dt_fac_max',
              'newton_max', 'it_target', 'use_it_target', 'it_p', 'cfl')
REQUIRED_MODEL_KEYS = ('d', 'rho_l0', 'ps', 'mul', 'mug', 'gas_eos', 'closure_riser',
                       'junction_model', 'dissipation_mode', 'flux_scheme')
CHANNELS = {
    'mdot_l_out_kg_s': ('liquid_outlet_mass_flow', 'kg/s'),
    'mdot_g_out_kg_s': ('gas_outlet_mass_flow', 'kg/s'),
    'mdot_l_base_kg_s': ('liquid_base_mass_flow', 'kg/s'),
    'mdot_g_base_kg_s': ('gas_base_mass_flow', 'kg/s'),
    'jl_out_m_s': ('liquid_outlet_superficial_velocity', 'm/s'),
    'jg_out_m_s': ('gas_outlet_superficial_velocity', 'm/s'),
    'Ml_kg': ('total_liquid_inventory', 'kg'),
    'Mg_kg': ('total_gas_inventory', 'kg'),
    'Ml_riser_kg': ('riser_liquid_inventory', 'kg'),
    'Mg_pipe_kg': ('pipeline_gas_inventory', 'kg'),
    'Mg_riser_kg': ('riser_gas_inventory', 'kg'),
    'Ml_pipe_kg': ('pipeline_liquid_inventory', 'kg'),
    'alpha_riser_base': ('riser_base_cell_void_fraction', '1'),
    'alpha_pipe_trace': ('pipeline_junction_trace_void_fraction', '1'),
    'alpha_riser_trace': ('riser_junction_trace_void_fraction', '1'),
    'dt_s': ('observed_accepted_time_step', 's'),
}


def same_number(a, b):
    x, y = finite(a), finite(b)
    if x is not None and y is not None:
        return math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-14)
    return a == b


def configuration_differences(a, b):
    pars_a, pars_b = a['parameters'], b['parameters']
    out = [key for key in MATCH_KEYS if not same_number(pars_a.get(key), pars_b.get(key))]
    if a['case'] != b['case']:
        out.append('case')
    return out


def discretization(record):
    pars = record['parameters']
    return (finite(pars.get('nodes')), finite(pars.get('n_riser')), finite(pars.get('dt_max_s')))


def source_signature(record):
    sources = record['manifest'].get('source_sha256')
    if not sources:
        return None
    return hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()


def mass_evidence(record):
    data = record.get('data')
    if not data:
        return record['upstream_metrics'].get('mass', {})
    result = {}
    for short, phase in [('l', 'liquid'), ('g', 'gas')]:
        error = [abs(x) for x in data.get('balance_' + short + '_kg', []) if x is not None]
        inventory = [x for x in data.get('M' + short + '_kg', []) if x is not None]
        influx = [x for x in data.get('M' + short + '_in_integral_kg', []) if x is not None]
        if error and inventory and influx:
            normalization = max(abs(inventory[0]), abs(influx[-1]), 1e-12)
            result[phase] = {'max_absolute_error_kg': max(error), 'normalization_kg': normalization,
                             'max_relative_error': max(error) / normalization,
                             'passes_1e_minus_6_relative': max(error) <= 1e-8 + 1e-6 * normalization}
    return result


def summarize_window(record, start, end, amplitude_floor):
    data = record.get('data')
    measures = {}
    warnings = []
    pkey = pressure_key(data, start, end)
    if pkey:
        measures['pressure'] = time_statistics(data, pkey, start, end) | {'unit': 'Pa', 'column': pkey}
    else:
        old = record['upstream_metrics']
        aligned = (same_number(old.get('window_start_s'), start) and same_number(old.get('window_end_s'), end))
        if aligned and old.get('cycles', {}).get('pressure_mean_Pa') is not None:
            cyc = old['cycles']
            pkey = old.get('pressure_location')
            measures['pressure'] = {'mean': finite(cyc.get('pressure_mean_Pa')),
                                    'std': finite(cyc.get('pressure_rms_Pa')),
                                    'peak_to_peak': finite(cyc.get('pressure_peak_to_peak_Pa')),
                                    'unit': 'Pa', 'column': pkey,
                                    'method': 'upstream metric; raw window unavailable; weighting not re-evaluated',
                                    'coverage_fraction': None}
            warnings.append('pressure_from_upstream_summary_only')
        else:
            warnings.append('pressure_unavailable_on_requested_window')
    for key, (name, unit) in CHANNELS.items():
        stats = time_statistics(data, key, start, end)
        if stats:
            measures[name] = stats | {'unit': unit, 'column': key}
            if stats['coverage_fraction'] < 1 - 1e-8:
                warnings.append(f'incomplete_channel_coverage:{key}')
    cycles = pressure_cycles(record, pkey, start, end, amplitude_floor)
    if cycles['median_period_s'] is None:
        warnings.append('no_resolved_period_on_requested_window')
    elif not cycles['ten_intervals_observed']:
        warnings.append('fewer_than_ten_pressure_intervals')
    if cycles.get('samples_per_period') is not None and cycles['samples_per_period'] < 20:
        warnings.append('fewer_than_twenty_observations_per_period')
    if cycles.get('period_cv') is not None and cycles['period_cv'] >= .2:
        warnings.append('periods_variable_or_modulated')
    blocks = []
    if pkey and data:
        for i in range(3):
            left, right = start + (end-start)*i/3, start + (end-start)*(i+1)/3
            stats = time_statistics(data, pkey, left, right)
            if stats:
                blocks.append({'start_s': left, 'end_s': right, 'mean_Pa': stats['mean'],
                               'std_Pa': stats['std'], 'peak_to_peak_Pa': stats['peak_to_peak']})
    if len(blocks) == 3:
        deviations = [block['std_Pa'] for block in blocks]
        if max(deviations) >= 1 and (max(deviations)-min(deviations))/max(deviations) > .20:
            warnings.append('pressure_std_differs_over20percent_among_window_thirds')
    return {'window_start_s': start, 'window_end_s': end, 'pressure_location': pkey,
            'measures': measures, 'cycles': cycles, 'pressure_window_thirds': blocks, 'warnings': warnings}


def make_pairs(records):
    pairs, excluded = [], []
    for left, right in itertools.combinations(records, 2):
        differences = configuration_differences(left, right)
        if differences:
            excluded.append({'runs': [left['run'], right['run']], 'reason': 'different_study_configuration',
                             'different_keys': differences})
            continue
        nl, rl, dl = discretization(left)
        nr, rr, dr = discretization(right)
        if any(value is None or value <= 0 for value in (nl, rl, dl, nr, rr, dr)):
            excluded.append({'runs': [left['run'], right['run']], 'reason': 'missing_discretization_metadata'})
            continue
        if same_number(dl, dr) and (nl != nr or rl != rr) and (nl - nr) * (rl - rr) >= 0:
            coarse, fine = (left, right) if nl + rl < nr + rr else (right, left)
            pairs.append(('mesh_at_same_dt_ceiling', coarse, fine))
        elif nl == nr and rl == rr and not same_number(dl, dr):
            coarse, fine = (left, right) if dl > dr else (right, left)
            pairs.append(('dt_ceiling_at_same_mesh', coarse, fine))
        else:
            excluded.append({'runs': [left['run'], right['run']],
                             'reason': 'both_mesh_and_dt_changed_or_no_refinement'})
    return pairs, excluded


def metric_differences(a, b):
    rows = []
    for measure in sorted(set(a['measures']) | set(b['measures'])):
        if measure == 'observed_accepted_time_step':
            statistics_keys = ('mean', 'minimum', 'maximum')
        elif 'inventory' in measure:
            statistics_keys = ('mean', 'std', 'peak_to_peak', 'change')
        elif measure == 'pressure':
            statistics_keys = ('mean', 'std', 'peak_to_peak')
        else:
            statistics_keys = ('mean', 'std', 'minimum', 'maximum')
        for statistic in statistics_keys:
            left, right = a['measures'].get(measure, {}), b['measures'].get(measure, {})
            va, vb = finite(left.get(statistic)), finite(right.get(statistic))
            comparable = (va is not None and vb is not None and
                          (measure != 'pressure' or a['pressure_location'] == b['pressure_location']) and
                          all(s.get('coverage_fraction') is None or s['coverage_fraction'] >= 1-1e-8 for s in (left, right)))
            delta = vb - va if comparable else None
            denominator = max(abs(va), abs(vb)) if delta is not None else None
            rows.append({'measure': measure, 'statistic': statistic, 'unit': left.get('unit', right.get('unit')),
                         'coarse': va, 'fine': vb, 'direct_difference_comparable': comparable, 'fine_minus_coarse': delta,
                         'absolute_difference': abs(delta) if delta is not None else None,
                         'relative_to_fine': abs(delta) / abs(vb) if delta is not None and vb != 0 else None,
                         'symmetric_relative_difference': abs(delta) / denominator if denominator else None})
    va, vb = a['cycles']['median_period_s'], b['cycles']['median_period_s']
    comparable = (va is not None and vb is not None and a['cycles']['method'] == b['cycles']['method'] and
                  a['pressure_location'] == b['pressure_location'])
    delta = vb - va if comparable else None
    rows.append({'measure': 'pressure_period', 'statistic': 'median', 'unit': 's', 'coarse': va, 'fine': vb,
                 'direct_difference_comparable': comparable, 'fine_minus_coarse': delta,
                 'absolute_difference': abs(delta) if delta is not None else None,
                 'relative_to_fine': abs(delta) / abs(vb) if delta is not None and vb else None,
                 'symmetric_relative_difference': abs(delta) / max(abs(va), abs(vb)) if delta is not None else None})
    return rows


def compare(records, start=900.0, end=None, amplitude_floor=100.0):
    runs = []
    for record in records:
        window = analysis_window(record, start, end)
        summary = summarize_window(record, *window, amplitude_floor) if window else None
        runs.append(json_ready(record) | {'analysis': summary, 'mass': mass_evidence(record),
                    'source_signature': source_signature(record),
                    'missing_model_metadata': [key for key in REQUIRED_MODEL_KEYS if key not in record['parameters']]})
    comparisons = []
    pairs, excluded = make_pairs(records)
    for kind, coarse, fine in pairs:
        flags = []
        if not coarse['completed'] or not fine['completed']:
            flags.append('incomplete_run_in_pair')
        signatures = (source_signature(coarse), source_signature(fine))
        if None in signatures:
            flags.append('implementation_provenance_missing')
        elif signatures[0] != signatures[1]:
            flags.append('different_implementation_sources')
        hashes = (coarse['manifest'].get('binary_sha256'), fine['manifest'].get('binary_sha256'))
        if all(hashes) and hashes[0] != hashes[1]:
            flags.append('different_executable_binaries')
        for record in (coarse, fine):
            if any(key not in record['parameters'] for key in REQUIRED_MODEL_KEYS):
                flags.append('physical_or_model_metadata_incomplete')
            mass = mass_evidence(record)
            if len(mass) != 2:
                flags.append('mass_balance_evidence_missing')
            elif not all(phase.get('passes_1e_minus_6_relative') for phase in mass.values()):
                flags.append('mass_balance_tolerance_failed')
        wa, wb = analysis_window(coarse, start, end), analysis_window(fine, start, end)
        item = {'kind': kind, 'coarse_run': coarse['run'], 'fine_run': fine['run'], 'case': coarse['case'],
                'coarse_discretization': dict(zip(('n_pipe', 'n_riser', 'dt_max_s'), discretization(coarse))),
                'fine_discretization': dict(zip(('n_pipe', 'n_riser', 'dt_max_s'), discretization(fine))),
                'automatic_convergence_claim': False, 'automatic_experimental_validation_claim': False}
        if not wa or not wb or min(wa[1], wb[1]) <= max(wa[0], wb[0]):
            flags.append('no_common_analysis_window')
            item.update(common_window_s=None, differences=[], status='insufficient_observation_window')
        else:
            common = (max(wa[0], wb[0]), min(wa[1], wb[1]))
            a = summarize_window(coarse, *common, amplitude_floor)
            b = summarize_window(fine, *common, amplitude_floor)
            flags += ['coarse:' + flag for flag in a['warnings']]
            flags += ['fine:' + flag for flag in b['warnings']]
            if a['pressure_location'] != b['pressure_location']:
                flags.append('pressure_locations_differ')
            if a['cycles']['method'] != b['cycles']['method']:
                flags.append('period_detection_methods_differ')
            pa, pb = a['measures'].get('pressure', {}), b['measures'].get('pressure', {})
            dynamic_status = 'insufficient_pressure_observation'
            if all(finite(p.get(key)) is not None for p in (pa, pb) for key in ('std', 'peak_to_peak')):
                resolved_a = pa['peak_to_peak'] >= amplitude_floor and pa['std'] >= 1
                resolved_b = pb['peak_to_peak'] >= amplitude_floor and pb['std'] >= 1
                changed = False
                if resolved_a != resolved_b:
                    flags.append('resolved_pressure_oscillation_presence_changes_under_refinement')
                    changed = True
                for key, minimum in [('std', 1.0), ('peak_to_peak', amplitude_floor)]:
                    scale = max(abs(pa[key]), abs(pb[key]))
                    if scale >= minimum and abs(pa[key]-pb[key])/scale > .10:
                        flags.append(f'pressure_{key}_differs_over10percent_under_refinement')
                        changed = True
                periods = (a['cycles']['median_period_s'], b['cycles']['median_period_s'])
                if all(value is not None for value in periods) and abs(periods[0]-periods[1])/max(periods) > .05:
                    flags.append('pressure_period_differs_over5percent_under_refinement')
                    changed = True
                dynamic_status = ('not_stable_across_refinement_controls' if changed else
                                  'both_pressure_signals_below_detection_floor' if not resolved_a and not resolved_b else
                                  'differences_within_declared_screening_thresholds; convergence_not_certified')
            if kind == 'dt_ceiling_at_same_mesh':
                actual_a = a['measures'].get('observed_accepted_time_step', {}).get('mean')
                actual_b = b['measures'].get('observed_accepted_time_step', {}).get('mean')
                if actual_a and actual_b and actual_b / actual_a > .9:
                    flags.append('dt_ceiling_changed_but_observed_mean_step_not_materially_reduced')
            item.update(common_window_s=list(common), coarse_statistics=a, fine_statistics=b,
                        differences=metric_differences(a, b), status='descriptive_comparison_available',
                        dynamic_comparison_status=dynamic_status)
        item['flags'] = sorted(set(flags))
        item['interpretation'] = ('Differences on a common window; no convergence order or GCI inferred. '
                                  'Mass conservation is necessary evidence and is not a convergence criterion by itself. '
                                  'Missing cycles prevent a period comparison, not a claim of physical stability.')
        comparisons.append(item)
    return {'schema': 'PRIMS-F.refinement-comparison.v1', 'created_utc': datetime.now(timezone.utc).isoformat(),
            'requested_window': {'start_s': start, 'end_s': end}, 'amplitude_floor_Pa': amplitude_floor,
            'statistical_method': 'Time-weighted piecewise-linear monitor integration on matched physical-time windows.',
            'time_control_warning': 'dt_max_s is an adaptive ceiling. Observed monitor dt samples do not contain every accepted step.',
            'convergence_policy': 'No automatic convergence or experimental validation claim; at least independent studies, adequate windows and physical evidence must be reviewed.',
            'dynamic_screening_policy': {'std_floor_Pa': 1.0, 'amplitude_floor_Pa': amplitude_floor,
                                         'amplitude_or_std_relative_difference': .10,
                                         'period_relative_difference': .05,
                                         'window_thirds_std_relative_spread': .20,
                                         'experimental_acceptance_thresholds': False},
            'runs': runs, 'comparisons': comparisons, 'excluded_pairs': excluded}


def markdown(report):
    lines = ['# Comparação de malha e passo temporal', '',
             '**Estas tabelas apresentam diferenças numéricas; não aprovam convergência ou captura experimental de slugging.**', '',
             f"Janela solicitada: {fmt(report['requested_window']['start_s'])} s até "
             f"{fmt(report['requested_window']['end_s']) if report['requested_window']['end_s'] is not None else 'o fim observado'}.", '',
             'As médias e desvios usam ponderação no tempo. `dt_max_s` é o teto de um passo adaptativo. '
             'A ausência de período resolvido não equivale a estabilidade física.', '',
             '| Execução | Concluída | Tempo alcançado (s) | Malha pipe+riser | dt máximo (s) | Pressão média (Pa) | ΔP (Pa) | Período (s) | Intervalos |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for run in report['runs']:
        pars, analysis = run['parameters'], run['analysis'] or {}
        pressure = analysis.get('measures', {}).get('pressure', {})
        cycles = analysis.get('cycles', {})
        lines.append('| ' + ' | '.join(map(fmt, [run['run'], run['completed'], run['reached_time_s'],
                     f"{fmt(pars.get('nodes'))}+{fmt(pars.get('n_riser'))}", pars.get('dt_max_s'),
                     pressure.get('mean'), pressure.get('peak_to_peak'), cycles.get('median_period_s'),
                     cycles.get('complete_intervals_found')])) + ' |')
    if not report['runs']:
        lines += ['', 'Nenhuma execução foi encontrada com o filtro solicitado.']
    if not report['comparisons']:
        lines += ['', 'Nenhum par independente de refinamento foi identificado com a mesma configuração física e numérica restante.']
    for pair in report['comparisons']:
        lines += ['', f"## {pair['kind']}: {pair['coarse_run']} → {pair['fine_run']}", '',
                  f"Janela comum: {pair['common_window_s']}; estado: `{pair['status']}`.", '',
                  f"Dinâmica: `{pair.get('dynamic_comparison_status', 'janela insuficiente')}`.", '',
                  'Limitações: ' + (', '.join(pair['flags']) if pair['flags'] else 'nenhuma das limitações automáticas listadas foi acionada; revisão continua necessária.'), '',
                  '| Medida | Estatística | Unidade | Menor resolução | Maior resolução | Diferença absoluta | Relativa à maior resolução |',
                  '|---|---|---|---:|---:|---:|---:|']
        for row in pair['differences']:
            relative = f"{100*row['relative_to_fine']:.5g}%" if row['relative_to_fine'] is not None else '—'
            lines.append('| ' + ' | '.join(map(fmt, [row['measure'], row['statistic'], row['unit'], row['coarse'],
                         row['fine'], row['absolute_difference'], relative])) + ' |')
    lines += ['', '## Limites da conclusão', '',
              '- Conservação de massa é verificada separadamente e não substitui concordância de pressão, amplitude, produção, inventários e período.',
              '- Duas resoluções fornecem uma diferença observada; este script não estima ordem de convergência nem GCI.',
              '- Execuções incompletas, configurações diferentes, canais ausentes, janelas incompatíveis e poucos ciclos permanecem sinalizados no JSON.',
              '- Os eventos físicos e seus limites de interpretação são tratados separadamente por `analyze_phase_sequence.py`.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metrics', type=Path)
    parser.add_argument('--results', type=Path, default=ROOT / 'validation/results')
    parser.add_argument('--run-dir', type=Path, action='append', default=[])
    parser.add_argument('--pattern', default='final_*')
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/analysis_final')
    parser.add_argument('--window-start', type=float, default=900.0)
    parser.add_argument('--window-end', type=float)
    parser.add_argument('--amplitude-floor', type=float, default=100.0)
    args = parser.parse_args()
    if args.window_start < 0 or args.amplitude_floor <= 0 or (args.window_end is not None and args.window_end <= args.window_start):
        parser.error('nonnegative start, positive detection floor and end>start are required')
    candidate = ROOT / 'validation/analysis_final/metrics.json'
    metrics = args.metrics or (candidate if candidate.is_file() else None)
    if args.metrics and not args.metrics.is_file():
        parser.error(f'metrics file not found: {args.metrics}')
    records = load_runs(args.results, args.pattern, metrics, args.run_dir)
    report = compare(records, args.window_start, args.window_end, args.amplitude_floor)
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / 'comparison.json', report)
    (args.output / 'comparison.md').write_text(markdown(report))
    print(json.dumps({'runs': len(records), 'pairs': len(report['comparisons']),
                      'completed': sum(r['completed'] for r in records), 'output': str(args.output)}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
