"""Analysis safeguards; synthetic signals are software fixtures, not benchmarks."""
import copy
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from validation_evidence import pressure_cycles, time_statistics
from compare_validation import compare, make_pairs
from analyze_phase_sequence import analyze, gas_restarts, inlet_reference


def fixture(name='flat', wave=False, cells=10, dt=.5, completed=True):
    area = math.pi * .05 ** 2 / 4
    times = [i / 5 for i in range(401)]
    data = {key: [] for key in ['p_junction_Pa', 'Ml_riser_kg', 'Mg_pipe_kg',
                                'jl_out_m_s', 'mdot_l_out_kg_s', 'mdot_g_out_kg_s', 'mdot_g_base_kg_s']}
    for t in times:
        u = t % 20
        if not wave:
            p, ml, mg, q, gb = 100000., 10., .1, .001, .001
        elif u < 5:
            p, ml, mg, q, gb = 100000 + 2000*u, 10 + .2*u, .1 + .001*u, .0005, 0.
        elif u < 10:
            p, ml, mg, q, gb = 110000., 11 - .05*(u-5), .105 + .0005*(u-5), .0035, 0.
        elif u < 15:
            p, ml, mg, q, gb = 110000 - 2000*(u-10), 10.75 - .15*(u-10), .1075 - .0015*(u-10), .001, .003
        else:
            p, ml, mg, q, gb = 100000., 10., .1, .0005, .001
        for key, value in [('p_junction_Pa', p), ('Ml_riser_kg', ml), ('Mg_pipe_kg', mg),
                           ('jl_out_m_s', q/area), ('mdot_l_out_kg_s', q*1000),
                           ('mdot_g_out_kg_s', gb), ('mdot_g_base_kg_s', gb)]:
            data[key].append(value)
    data['time_s'] = times
    data['dt_s'] = [dt] * len(times)
    for key, value in [('Ml_kg', 50.), ('Mg_kg', .2), ('balance_l_kg', 0.), ('balance_g_kg', 0.)]:
        data[key] = [value] * len(times)
    data['Ml_in_integral_kg'] = times[:]
    data['Mg_in_integral_kg'] = [t*.001 for t in times]
    return {'run': name, 'case': 'SYNTHETIC_SOFTWARE_FIXTURE', 'folder': '/synthetic/not-a-benchmark',
            'parameters': {'nodes': cells, 'n_riser': cells, 'dt_max_s': dt, 'd': .05,
                           'rho_l0': 1000., 'ps': 100000., 'jl': .001/area,
                           'mul': 1/area, 'mug': .001/area, 'gas_eos': 'BAROTROPIC',
                           'closure_riser': 3, 'junction_model': 'MATCHED',
                           'dissipation_mode': 'LAGGED', 'flux_scheme': 'ROE'},
            'manifest': {'source_sha256': {'synthetic': 'same'}, 'binary_sha256': 'same'},
            'completed': completed, 'termination_status': 0 if completed else 6000,
            'reached_time_s': 80., 'requested_time_s': 80., 'data': data,
            'upstream_metrics': {}, 'warnings': [], 'provenance': {}}


class EvidenceTests(unittest.TestCase):
    def test_time_weighted_statistics_on_irregular_samples(self):
        stats = time_statistics({'time_s': [0., 1., 10.], 'x': [0., 1., 10.]}, 'x', 0, 10)
        self.assertAlmostEqual(stats['mean'], 5)
        self.assertAlmostEqual(stats['std'], math.sqrt(100/12))
        self.assertAlmostEqual(stats['coverage_fraction'], 1)

    def test_missing_segments_are_not_zero_or_interpolated_away(self):
        stats = time_statistics({'time_s': [0., 1., 2., 3.], 'x': [0., 1., None, 3.]}, 'x', 0, 3)
        self.assertAlmostEqual(stats['covered_duration_s'], 1)
        self.assertAlmostEqual(stats['coverage_fraction'], 1/3)

    def test_only_independent_refinement_pairs(self):
        pairs, excluded = make_pairs([fixture('base'), fixture('mesh', cells=20), fixture('time', dt=.25)])
        self.assertEqual({p[0] for p in pairs}, {'mesh_at_same_dt_ceiling', 'dt_ceiling_at_same_mesh'})
        self.assertEqual(len(pairs), 2)
        self.assertEqual(len(excluded), 1)

    def test_different_physics_is_not_a_refinement_pair(self):
        a, b = fixture('base'), fixture('different', cells=20)
        b['parameters']['ps'] = 200000
        pairs, excluded = make_pairs([a, b])
        self.assertFalse(pairs)
        self.assertIn('ps', excluded[0]['different_keys'])

    def test_mass_pass_does_not_establish_convergence(self):
        report = compare([fixture('a'), fixture('b', cells=20, completed=False)], start=0)
        pair = report['comparisons'][0]
        self.assertFalse(pair['automatic_convergence_claim'])
        self.assertIn('incomplete_run_in_pair', pair['flags'])
        self.assertIn('coarse:no_resolved_period_on_requested_window', pair['flags'])

    def test_new_source_hash_is_flagged(self):
        a, b = fixture('a'), fixture('b', cells=20)
        b['manifest']['source_sha256']['synthetic'] = 'changed'
        pair = compare([a, b], start=0)['comparisons'][0]
        self.assertIn('different_implementation_sources', pair['flags'])

    def test_unreached_window_never_falls_back_to_a_pilot_window(self):
        report = compare([fixture('a'), fixture('b', cells=20)], start=900)
        self.assertIsNone(report['runs'][0]['analysis'])
        self.assertIn('no_common_analysis_window', report['comparisons'][0]['flags'])

    def test_flat_signals_create_no_periods_or_phase_sequences(self):
        record = fixture()
        cycles = pressure_cycles(record, 'p_junction_Pa', 0, 80)
        self.assertIsNone(cycles['median_period_s'])
        out = analyze(record, start=0)
        self.assertEqual(out['classification'], 'no_resolved_events_at_declared_thresholds')
        self.assertFalse(out['ordered_candidates'])
        self.assertTrue(all(not events for events in out['events'].values()))

    def test_missing_compartment_inventory_blocks_ordered_claim(self):
        record = fixture(wave=True)
        del record['data']['Mg_pipe_kg']
        out = analyze(record, start=0)
        self.assertEqual(out['classification'], 'insufficient_named_signals')
        self.assertIn('gas_pipeline_inventory', out['missing_signals'])
        self.assertFalse(out['ordered_candidates'])

    def test_synthetic_order_is_detected_but_not_validated(self):
        out = analyze(fixture(wave=True), start=0)
        self.assertGreaterEqual(len(out['ordered_candidates']), 3)
        self.assertFalse(out['automatic_slugging_validation'])
        self.assertFalse(out['automatic_flow_regime_classification'])
        self.assertTrue(out['repeated_ordered_candidates'])

    def test_monotone_pressure_is_not_a_cycle(self):
        record = fixture()
        record['data']['p_junction_Pa'] = [100000+2000*t for t in record['data']['time_s']]
        cycles = pressure_cycles(record, 'p_junction_Pa', 0, 80)
        self.assertIsNone(cycles['median_period_s'])
        self.assertEqual(cycles['complete_intervals_found'], 0)

    def test_separate_short_low_flux_visits_are_not_a_continuous_dwell(self):
        events = gas_restarts([0, 1, 2, 3, 4], [0, .15, .15, 0, .5], .05, .25, .5)
        self.assertFalse(events)

    def test_effective_mass_override_precedes_nominal_volumetric_input(self):
        record = fixture()
        record['parameters'].update(d=1., jl=1., mul=2000., rho_l0=1000.)
        value, definition = inlet_reference(record, record['data'], 0, 80, 'l', True)
        self.assertAlmostEqual(value, math.pi/2)
        self.assertIn('imposed_liquid_mass', definition)

    def test_different_numeric_controller_is_not_a_one_factor_study(self):
        a, b = fixture('a'), fixture('b', cells=20)
        a['parameters']['dt_growth'] = 1.2
        b['parameters']['dt_growth'] = 1.5
        pairs, excluded = make_pairs([a, b])
        self.assertFalse(pairs)
        self.assertIn('dt_growth', excluded[0]['different_keys'])

    def test_amplitude_change_is_flagged_even_if_mass_balances_pass(self):
        pair = compare([fixture('flat'), fixture('oscillating', wave=True, dt=.125)], start=0)['comparisons'][0]
        self.assertIn('resolved_pressure_oscillation_presence_changes_under_refinement', pair['flags'])
        self.assertEqual(pair['dynamic_comparison_status'], 'not_stable_across_refinement_controls')
        self.assertFalse(pair['automatic_convergence_claim'])

    def test_pressure_oscillation_with_positive_gas_is_preserved_as_positive_evidence(self):
        record = fixture(wave=True)
        record['data']['mdot_g_base_kg_s'] = [.002] * len(record['data']['time_s'])
        out = analyze(record, start=0)
        self.assertEqual(out['classification'], 'pressure_oscillations_with_observed_positive_base_gas_flow')
        self.assertFalse(out['ordered_candidates'])
        self.assertTrue(out['gas_base_observations']['no_restart_does_not_exclude_osc_or_ss3'])


if __name__ == '__main__':
    unittest.main()
