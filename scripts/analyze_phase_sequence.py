#!/usr/bin/env python3
"""Measure qualitative inventory/production/gas-transfer sequences in named signals.

These are numerical signal proxies, not observations of a moving slug front.
Flat signals do not generate events. Missing compartment inventories or base
gas flow prevent an ordered-sequence claim. Absence of this particular proxy
does not exclude SS3, OSC or hydrodynamic slug flow.
"""
from __future__ import annotations

import argparse
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

from validation_evidence import (ROOT, analysis_window, clipped_series, finite, fmt, json_ready,
                                 load_runs, median_interval, pressure_cycles, pressure_key,
                                 time_statistics, value_at, write_json)

DEFAULTS = {'trend_window_s': 2.0, 'min_event_duration_s': .5, 'max_sequence_span_s': 300.0,
            'pressure_change_floor_Pa': 100.0, 'inventory_fraction_of_range': .005,
            'liquid_inventory_floor_kg': 1e-5, 'gas_inventory_floor_kg': 1e-7,
            'liquid_production_floor_m3_s': 1e-8, 'liquid_production_floor_kg_s': 1e-5,
            'production_excess_fraction_of_reference': .10,
            'production_excursion_fraction_of_reference': .05,
            'gas_base_floor_kg_s': 1e-9, 'gas_low_fraction_of_reference': .05,
            'gas_high_fraction_of_reference': .25, 'pressure_cycle_amplitude_floor_Pa': 100.0}


def episodes(times, mask, name, minimum_duration):
    result = []
    first = None
    for i in range(len(times) + 1):
        active = i < len(times) and bool(mask[i])
        if active and first is None:
            first = i
        if not active and first is not None:
            last = i - 1
            if times[last] - times[first] >= minimum_duration:
                result.append({'id': f'{name}_{len(result)+1}', 'start_s': times[first], 'end_s': times[last],
                               'duration_s': times[last] - times[first], 'observations': last - first + 1,
                               'evidence_type': 'thresholded_numerical_signal_episode'})
            first = None
    return result


def trend_changes(data, key, times, start, end, span):
    result = []
    for time in times:
        left, right = max(start, time - span / 2), min(end, time + span / 2)
        a, b = value_at(data['time_s'], data[key], left), value_at(data['time_s'], data[key], right)
        result.append(b - a if a is not None and b is not None and right > left else None)
    return result


def change(data, key, event):
    stats = time_statistics(data, key, event['start_s'], event['end_s'])
    return stats['change'] if stats else None


def inlet_reference(record, data, start, end, phase, volumetric=False):
    pars = record['parameters']
    diameter = finite(pars.get('d'))
    area = math.pi * diameter ** 2 / 4 if diameter is not None and diameter > 0 else None
    mass_flux = finite(pars.get('mu' + phase))
    if area and mass_flux is not None and mass_flux >= 0:
        mass_flow = area * mass_flux
        if not volumetric:
            return mass_flow, 'A*imposed_phase_mass_flux from metadata'
        density = finite(pars.get('rho_l0'))
        if density and density > 0:
            return mass_flow / density, 'imposed_liquid_mass_flow/rho_l0; reference-volume proxy'
    if volumetric and phase == 'l' and area and finite(pars.get('jl')) is not None:
        return area * finite(pars['jl']), 'A*jl_reference from metadata; no effective mass override available'
    integral = 'M' + phase + '_in_integral_kg'
    if integral in data:
        a, b = value_at(data['time_s'], data[integral], start), value_at(data['time_s'], data[integral], end)
        if a is not None and b is not None and end > start and not volumetric:
            return (b - a) / (end - start), 'mean influx from cumulative inlet mass over analysis window'
    return None, 'reference unavailable'


def liquid_production_series(record, data):
    if 'Ql_out_m3_s' in data:
        return 'Ql_out_m3_s', 'm3/s', 'named volumetric outlet flow', True
    diameter = finite(record['parameters'].get('d'))
    if 'jl_out_m_s' in data and diameter and diameter > 0:
        area = math.pi * diameter ** 2 / 4
        data['derived_Ql_out_m3_s'] = [area * value if value is not None else None for value in data['jl_out_m_s']]
        return 'derived_Ql_out_m3_s', 'm3/s', 'A*jl_out, circular area from metadata D', True
    if 'mdot_l_out_kg_s' in data:
        return 'mdot_l_out_kg_s', 'kg/s', 'mass-production proxy; volumetric outlet flow unavailable', False
    return None, None, 'liquid production signal unavailable', False


def gas_restarts(times, values, low, high, minimum_duration):
    events = []
    first_low = last_low = None
    reversed_flow = False
    qualified = None
    for i, value in enumerate(values):
        if value is None:
            first_low = last_low = None
            reversed_flow = False
            qualified = None
            continue
        if value <= low:
            if first_low is None:
                first_low = times[i]
            last_low = times[i]
            reversed_flow |= value < 0
            if last_low - first_low >= minimum_duration:
                qualified = (first_low, last_low, reversed_flow)
        else:
            # The dwell interval must be contiguous below the LOW threshold.
            # A completed dwell may arm the later hysteresis upcrossing, but
            # separated short visits are never summed into a fictitious dwell.
            first_low = last_low = None
            reversed_flow = False
        if value >= high:
            if qualified is not None:
                observed_start, observed_end, observed_reverse = qualified
                crossing = times[i]
                if i and values[i - 1] is not None and values[i - 1] < high and value > values[i - 1]:
                    crossing = times[i - 1] + (times[i] - times[i - 1]) * (high - values[i - 1]) / (value - values[i - 1])
                events.append({'id': f'gas_restart_{len(events)+1}', 'time_s': crossing,
                               'preceding_low_flux_start_s': observed_start, 'last_low_flux_time_s': observed_end,
                               'observed_low_flux_duration_s': observed_end - observed_start,
                               'preceding_reverse_flow_observed': observed_reverse,
                               'low_threshold_kg_s': low, 'high_threshold_kg_s': high,
                               'evidence_type': 'base_gas_flux_restart_proxy; not a measured front position'})
            first_low = last_low = None
            reversed_flow = False
            qualified = None
    return events


def ordered_sequences(accumulations, productions, restarts, blowdowns, span, overlap):
    chains, used_restarts = [], set()
    for accumulation in accumulations:
        start = accumulation['start_s']
        for production in productions:
            if production['start_s'] < start - overlap or production['end_s'] < accumulation['end_s'] - overlap:
                continue
            if production['start_s'] > start + span:
                break
            matched = None
            for restart in restarts:
                if restart['id'] in used_restarts:
                    continue
                if restart['time_s'] < max(accumulation['end_s'], production['start_s']) - overlap:
                    continue
                if restart['time_s'] > start + span:
                    break
                blowdown = next((event for event in blowdowns
                                 if event['end_s'] >= restart['time_s'] and
                                 event['start_s'] >= restart['time_s'] - overlap and
                                 event['end_s'] <= start + span), None)
                if blowdown:
                    matched = (restart, blowdown)
                    break
            if matched:
                restart, blowdown = matched
                used_restarts.add(restart['id'])
                chains.append({'id': f'ordered_episode_{len(chains)+1}',
                               'start_s': start, 'end_s': blowdown['end_s'],
                               'accumulation_id': accumulation['id'], 'liquid_production_id': production['id'],
                               'gas_restart_id': restart['id'], 'blowdown_id': blowdown['id'],
                               'allowed_overlap_s': overlap,
                               'evidence_type': 'ordered numerical proxy; not an experimentally confirmed slug cycle'})
                break
    return chains


def analyze(record, start=900.0, end=None, thresholds=None):
    settings = dict(DEFAULTS)
    settings.update(thresholds or {})
    output = json_ready(record)
    output['warnings'] = list(record.get('warnings', []))
    output.update(events={}, ordered_candidates=[], automatic_slugging_validation=False,
                  automatic_flow_regime_classification=False, thresholds=settings,
                  ordered_sequence_scope=('Strong-blockage/restart motif, useful for an SS1-like narrative; '
                                          'not a necessary condition for OSC or SS3, which may retain gas penetration.'),
                  interpretation=('Inventories, pressure and phase fluxes are numerical observables. '
                                  'Compression, penetration and blowdown labels are qualitative proxies. '
                                  'No moving-front position or experimental phase timing was measured. '
                                  'Absence of this ordered proxy does not exclude SS3, OSC or hydrodynamic slug flow.'))
    window = analysis_window(record, start, end)
    if not window or not record.get('data'):
        output.update(classification='insufficient_observation_window', window=None,
                      missing_signals=['raw_monitor' if not record.get('data') else 'requested_time_window'])
        return output
    start, end = window
    data = dict(record['data'])
    pkey = pressure_key(data, start, end)
    qkey, qunit, qdefinition, volumetric = liquid_production_series(record, data)
    channels = {'pressure': pkey, 'liquid_riser_inventory': 'Ml_riser_kg',
                'gas_pipeline_inventory': 'Mg_pipe_kg', 'liquid_production': qkey,
                'gas_base_mass_flow': 'mdot_g_base_kg_s'}
    measures, missing = {}, []
    for label, key in channels.items():
        stats = time_statistics(data, key, start, end) if key else None
        if not stats or stats['coverage_fraction'] < 1 - 1e-8:
            missing.append(label)
        if stats:
            measures[label] = stats | {'column': key}
    output.update(window={'start_s': start, 'end_s': end}, measured_signals=measures,
                  signal_columns=channels, liquid_production_unit=qunit,
                  liquid_production_definition=qdefinition, missing_signals=missing,
                  pressure_cycles=pressure_cycles(record, pkey, start, end, settings['pressure_cycle_amplitude_floor_Pa']))
    if missing:
        output['classification'] = 'insufficient_named_signals'
        return output
    times, _ = clipped_series(data, pkey, start, end)
    output['median_observation_interval_s'] = median_interval(times)
    if len(times) < 4:
        output['classification'] = 'insufficient_samples'
        return output
    pstats, mlstats, mgstats = (measures[key] for key in ('pressure', 'liquid_riser_inventory', 'gas_pipeline_inventory'))
    pchange = settings['pressure_change_floor_Pa']
    mlchange = max(settings['liquid_inventory_floor_kg'], settings['inventory_fraction_of_range'] * mlstats['peak_to_peak'])
    mgchange = max(settings['gas_inventory_floor_kg'], settings['inventory_fraction_of_range'] * mgstats['peak_to_peak'])
    span = settings['trend_window_s']
    dp = trend_changes(data, pkey, times, start, end, span)
    dml = trend_changes(data, 'Ml_riser_kg', times, start, end, span)
    dmg = trend_changes(data, 'Mg_pipe_kg', times, start, end, span)
    accumulation_mask = [all(value is not None for value in values) and
                         values[0] > pchange and values[1] > mlchange and values[2] > mgchange
                         for values in zip(dp, dml, dmg)]
    accumulations = episodes(times, accumulation_mask, 'accumulation', settings['min_event_duration_s'])
    for event in accumulations:
        event['pressure_rise_Pa'] = change(data, pkey, event)
        event['liquid_riser_gain_kg'] = change(data, 'Ml_riser_kg', event)
        event['gas_pipeline_gain_kg'] = change(data, 'Mg_pipe_kg', event)
        event['interpretation'] = 'Coincident accumulation and pressure rise; compression is an inferred proxy, not a front measurement.'
    qref, qref_definition = inlet_reference(record, data, start, end, 'l', volumetric)
    qstats = measures['liquid_production']
    if qref is None:
        _, window_values = clipped_series(data, qkey, start, end)
        qref = statistics.median(value for value in window_values if value is not None)
        qref_definition = 'median observed outlet signal; imposed inlet reference unavailable'
    qfloor = settings['liquid_production_floor_m3_s' if volumetric else 'liquid_production_floor_kg_s']
    qhigh = qref + max(qfloor, settings['production_excess_fraction_of_reference'] * abs(qref), .2 * qstats['peak_to_peak'])
    qrange_floor = max(qfloor, settings['production_excursion_fraction_of_reference'] * abs(qref))
    qvalues = [value_at(data['time_s'], data[qkey], time) for time in times]
    production_mask = [value is not None and value >= qhigh and qstats['peak_to_peak'] >= qrange_floor for value in qvalues]
    productions = episodes(times, production_mask, 'liquid_production', settings['min_event_duration_s'])
    for event in productions:
        stats = time_statistics(data, qkey, event['start_s'], event['end_s'])
        event.update(peak=stats['maximum'], mean=stats['mean'], unit=qunit,
                     threshold=qhigh, reference=qref, reference_definition=qref_definition)
    gref, gref_definition = inlet_reference(record, data, start, end, 'g')
    gas_values = [value_at(data['time_s'], data['mdot_g_base_kg_s'], time) for time in times]
    if gref is None or gref <= 0:
        positives = sorted(value for value in gas_values if value is not None and value > 0)
        gref = positives[int(.9 * (len(positives)-1))] if positives else 0.0
        gref_definition = '90th percentile of observed positive base gas flow; inlet reference unavailable'
    glow = max(settings['gas_base_floor_kg_s'], settings['gas_low_fraction_of_reference'] * gref)
    ghigh = max(2 * glow, settings['gas_high_fraction_of_reference'] * gref)
    restarts = gas_restarts(times, gas_values, glow, ghigh, settings['min_event_duration_s'])
    blowdown_mask = [pd is not None and gv is not None and pd < -pchange and gv >= ghigh
                     for pd, gv in zip(dp, gas_values)]
    blowdowns = episodes(times, blowdown_mask, 'blowdown', settings['min_event_duration_s'])
    for event in blowdowns:
        event['pressure_change_Pa'] = change(data, pkey, event)
        event['gas_pipeline_change_kg'] = change(data, 'Mg_pipe_kg', event)
        event['mean_gas_base_kg_s'] = time_statistics(data, 'mdot_g_base_kg_s', event['start_s'], event['end_s'])['mean']
        gout = time_statistics(data, 'mdot_g_out_kg_s', event['start_s'], event['end_s'])
        if gout:
            event['mean_gas_outlet_kg_s'] = gout['mean']
            event['gas_outlet_exceeds_reference_inlet'] = gout['mean'] > gref
        event['interpretation'] = 'Pressure fall with positive base gas transfer; blowdown candidate, not direct experimental confirmation.'
    chains = ordered_sequences(accumulations, productions, restarts, blowdowns,
                               settings['max_sequence_span_s'], settings['trend_window_s'])
    output.update(events={'accumulation_compression_proxies': accumulations,
                          'liquid_production_excursions': productions,
                          'base_gas_restart_proxies': restarts, 'blowdown_proxies': blowdowns},
                  ordered_candidates=chains,
                  resolved_thresholds={'pressure_change_Pa': pchange, 'riser_liquid_change_kg': mlchange,
                                       'pipeline_gas_change_kg': mgchange, 'production_high': qhigh,
                                       'production_reference': qref, 'production_reference_definition': qref_definition,
                                       'production_excursion_floor': qrange_floor, 'production_unit': qunit,
                                       'gas_reference_kg_s': gref, 'gas_reference_definition': gref_definition,
                                       'gas_low_kg_s': glow, 'gas_high_kg_s': ghigh},
                  flat_or_unresolved_signals={'pressure': pstats['peak_to_peak'] < pchange,
                                             'liquid_riser_inventory': mlstats['peak_to_peak'] < mlchange,
                                             'gas_pipeline_inventory': mgstats['peak_to_peak'] < mgchange,
                                             'liquid_production': qstats['peak_to_peak'] < qrange_floor},
                  repeated_ordered_candidates=len(chains) >= 2)
    output['gas_base_observations'] = {
        'minimum_kg_s': measures['gas_base_mass_flow']['minimum'],
        'strictly_positive_at_observed_times': measures['gas_base_mass_flow']['minimum'] > 0,
        'above_low_flow_threshold_through_observation': measures['gas_base_mass_flow']['minimum'] > glow,
        'low_flux_detector_definition': 'At or below5% of reference gas influx (subject to absolute floor), continuously observed for at least configured minimum duration; restart above25% uses hysteresis.',
        'no_restart_does_not_exclude_osc_or_ss3': True}
    output['pressure_window_thirds'] = []
    for i in range(3):
        left, right = start+(end-start)*i/3, start+(end-start)*(i+1)/3
        stats = time_statistics(data, pkey, left, right)
        output['pressure_window_thirds'].append({'start_s': left, 'end_s': right,
                                                 'std_Pa': stats['std'], 'peak_to_peak_Pa': stats['peak_to_peak']})
    if chains:
        output['classification'] = 'ordered_episode_candidates' if len(chains) > 1 else 'single_ordered_episode_candidate'
    elif any(output['events'].values()):
        output['classification'] = 'unmatched_or_partial_transient_events'
    else:
        output['classification'] = 'no_resolved_events_at_declared_thresholds'
    if not chains and output['pressure_cycles']['median_period_s'] is not None and measures['gas_base_mass_flow']['minimum'] > 0:
        output['classification'] = 'pressure_oscillations_with_observed_positive_base_gas_flow'
    if not record['completed']:
        output['warnings'].append('incomplete_simulation_limits_sequence_observation')
    if any(value is not None and value < 0 for value in gas_values):
        output['warnings'].append('reverse_gas_base_flow_observed; low-flux proxy is not necessarily blockage')
    barrier = time_statistics(data, 'blocking_active', start, end)
    if barrier:
        output['model_imposed_barrier_active_time_fraction'] = barrier['mean']
    if str(record['parameters'].get('junction_model', '')).upper() == 'BLOCKING':
        output['warnings'].append('BLOCKING is an imposed numerical/model barrier; detected ordering does not validate a moving-front model')
    return output


def markdown(report):
    lines = ['# Sequência física: medidas e indicadores qualitativos', '',
             '**Os eventos abaixo são indicadores de sinais numéricos. Não constituem comprovação experimental de slugging nem observação direta de uma frente de líquido.**', '',
             'A identificação combina aumento de líquido no riser, ganho de gás no pipeline e elevação de pressão; '
             'produção de líquido acima da referência; retomada do fluxo de gás na base após fluxo baixo; '
             'e queda de pressão com transferência positiva de gás. Os limiares são escolhas declaradas de análise.', '',
             'A cadeia com retomada de gás é um indicador de **bloqueio forte e liberação**, útil para uma narrativa do tipo SS1. '
             'O detector exige fluxo gasoso na base abaixo de 5% da entrada de referência durante pelo menos 0,5 s '
             '(valores padrão), seguido de subida acima de 25%. **OSC e SS3 podem manter passagem de gás:** '
             'zero cadeias desse tipo não prova ausência dessas instabilidades. Oscilações com gás observado sempre '
             'positivo são registradas como um resultado próprio.', '',
             '| Execução | Concluída | Classificação descritiva | Acumulação | Produção | Retomada de gás | Despressurização | Sequências ordenadas |',
             '|---|---|---|---:|---:|---:|---:|---:|']
    for run in report['runs']:
        events = run.get('events', {})
        values = [run['run'], run['completed'], run['classification'],
                  len(events.get('accumulation_compression_proxies', [])), len(events.get('liquid_production_excursions', [])),
                  len(events.get('base_gas_restart_proxies', [])), len(events.get('blowdown_proxies', [])),
                  len(run.get('ordered_candidates', []))]
        lines.append('| ' + ' | '.join(map(fmt, values)) + ' |')
    if not report['runs']:
        lines += ['', 'Nenhuma execução corresponde ao filtro solicitado.']
    for run in report['runs']:
        lines += ['', '## ' + run['run'], '', f"Janela: {run.get('window')}; tempo alcançado: {fmt(run['reached_time_s'])} s.", '',
                  'Canais: `' + str(run.get('signal_columns', {})) + '`.', '',
                  'Canais insuficientes: ' + (', '.join(run.get('missing_signals', [])) or 'nenhum canal obrigatório ausente.'), '',
                  'Limitações: ' + (', '.join(run.get('warnings', [])) or 'nenhuma limitação adicional de leitura; a interpretação permanece qualitativa.'), '']
        if run.get('resolved_thresholds'):
            lines += ['Limiares efetivos:', '', '| Grandeza | Valor |', '|---|---|']
            lines += [f'| {fmt(key)} | {fmt(value)} |' for key, value in run['resolved_thresholds'].items()]
        if run.get('pressure_window_thirds'):
            lines += ['', '| Terço da janela (s) | Desvio de pressão (Pa) | ΔP (Pa) |', '|---|---:|---:|']
            lines += [f"| {fmt(b['start_s'])}–{fmt(b['end_s'])} | {fmt(b['std_Pa'])} | {fmt(b['peak_to_peak_Pa'])} |"
                      for b in run['pressure_window_thirds']]
        if run.get('ordered_candidates'):
            lines += ['', '| Episódio candidato | Início (s) | Fim (s) | Acumulação | Produção | Retomada de gás | Despressurização |',
                      '|---|---:|---:|---|---|---|---|']
            for event in run['ordered_candidates']:
                lines.append('| ' + ' | '.join(map(fmt, [event[key] for key in ('id', 'start_s', 'end_s',
                             'accumulation_id', 'liquid_production_id', 'gas_restart_id', 'blowdown_id')])) + ' |')
    lines += ['', '## Alcance', '',
              '- Os valores de inventários, pressão e vazões provêm de canais nomeados do solver. Compressão, penetração e blowdown são interpretações indiretas dessas medidas.',
              '- Sinais planos não recebem ciclos artificiais. Fluxo de gás inicialmente positivo não é contado como uma nova penetração sem um intervalo anterior de fluxo baixo.',
              '- Uma excursão isolada é um transiente candidato; a tabela não transforma o número de episódios em ciclos experimentais comprovados.',
              '- Esta sequência é um indicador particular. Sua ausência não exclui SS3, OSC ou slug flow hidrodinâmico; mecanismos diferentes não devem ser forçados a uma sequência previamente escolhida.',
              '- Falta de dados espaciais, duração insuficiente, domínio da correlação e limites da junção continuam relevantes mesmo quando existe uma sequência ordenada.', '']
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
    parser.add_argument('--trend-window-s', type=float, default=DEFAULTS['trend_window_s'])
    parser.add_argument('--min-event-duration-s', type=float, default=DEFAULTS['min_event_duration_s'])
    parser.add_argument('--max-sequence-span-s', type=float, default=DEFAULTS['max_sequence_span_s'])
    parser.add_argument('--pressure-change-floor-pa', type=float, default=DEFAULTS['pressure_change_floor_Pa'])
    args = parser.parse_args()
    if args.window_start < 0 or (args.window_end is not None and args.window_end <= args.window_start):
        parser.error('window must have nonnegative start and end>start')
    settings = {'trend_window_s': args.trend_window_s, 'min_event_duration_s': args.min_event_duration_s,
                'max_sequence_span_s': args.max_sequence_span_s, 'pressure_change_floor_Pa': args.pressure_change_floor_pa}
    if min(settings.values()) <= 0:
        parser.error('positive event-analysis controls are required')
    candidate = ROOT / 'validation/analysis_final/metrics.json'
    metrics = args.metrics or (candidate if candidate.is_file() else None)
    if args.metrics and not args.metrics.is_file():
        parser.error(f'metrics file not found: {args.metrics}')
    records = load_runs(args.results, args.pattern, metrics, args.run_dir)
    report = {'schema': 'PRIMS-F.phase-sequence.v1', 'created_utc': datetime.now(timezone.utc).isoformat(),
              'requested_window': {'start_s': args.window_start, 'end_s': args.window_end},
              'thresholds_are_analysis_choices_not_experimental_calibration': True,
              'runs': [analyze(record, args.window_start, args.window_end, settings) for record in records]}
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / 'phase_sequence.json', report)
    (args.output / 'phase_sequence.md').write_text(markdown(report))
    print(f"runs={len(records)} output={args.output}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
