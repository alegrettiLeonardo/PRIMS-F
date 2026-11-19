import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('accepted_mass', ROOT / 'scripts/analyze_accepted_conservation.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def fixture():
    """Signed synthetic transport with a known, deliberately nonzero defect."""
    rows = []
    mass = {'l': [10.0, 20.0], 'g': [1.0, 2.0]}
    initial = {p: values[:] for p, values in mass.items()}
    accum = {(p, f): 0.0 for p in audit.PHASES for f in audit.FACES}
    abssum = accum.copy()
    residual_sum = {(p, r): 0.0 for p in audit.PHASES for r in audit.REGIONS}
    time = 0.0
    for step, dt in enumerate((0.0, 0.1, 0.2)):
        time += dt
        row = dict(step=step, time_s=time, dt_s=dt, residual_norm=1e-8 if step else 0, flux_status=0)
        for p in audit.PHASES:
            flows = [2., 1., .5, .5] if p == 'l' else [-.2, -.4, -.1, -.1]
            if step == 2:
                flows = [-value for value in flows]
            defects = [1e-4, -2e-4] if step == 1 else [0., 0.]
            old = mass[p][:]
            mass[p][0] += dt * (flows[0] - flows[2]) + defects[0]
            mass[p][1] += dt * (flows[3] - flows[1]) + defects[1]
            for f, value in zip(audit.FACES, flows):
                row[f'mdot_{p}_{f}_kg_s'] = value
                accum[p, f] += dt * value
                abssum[p, f] += dt * abs(value)
                row[f'M{p}_{f}_integral_kg'] = accum[p, f]
                row[f'M{p}_{f}_abs_integral_kg'] = abssum[p, f]
            for r, value, value0, defect in zip(audit.REGIONS, [sum(mass[p]), *mass[p]],
                                              [sum(initial[p]), *initial[p]], [sum(defects), *defects]):
                row[f'M{p}{r}_kg'] = value
                incoming, outgoing = {'': ('in', 'out'), '_pipe': ('in', 'base_left'), '_riser': ('base_right', 'out')}[r]
                row[f'balance_{p}{r}_kg'] = value - value0 - accum[p, incoming] + accum[p, outgoing]
                row[f'step_balance_{p}{r}_kg'] = defect
                row[f'mass_residual_{p}{r}_kg'] = defect
                residual_sum[p, r] += defect
                row[f'mass_residual_sum_{p}{r}_kg'] = residual_sum[p, r]
        rows.append(row)
    return rows


class AcceptedConservationTests(unittest.TestCase):
    def test_status_retains_textual_termination_reason(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'run_status.txt'
            path.write_text('status=6000\ntermination_reason=LINE_SEARCH_EXHAUSTED\naccepted_steps=3\n')
            result = audit.read_status(path)
            self.assertEqual(result['status'], 6000)
            self.assertEqual(result['termination_reason'], 'LINE_SEARCH_EXHAUSTED')
            self.assertEqual(result['accepted_steps'], 3)

    def analyze(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'conservation.csv'
            with path.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=fixture()[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            return audit.audit_ledger(path)

    def test_signed_transport_and_nonzero_newton_defect(self):
        result = self.analyze(fixture())
        self.assertTrue(result['arithmetic_contract_passes'], result['mismatches'])
        self.assertEqual(result['accepted_steps'], 2)
        self.assertFalse(result['balances']['l_system']['passes_inventory_1e6'])
        self.assertAlmostEqual(result['balances']['g_system']['absolute_processed_mass_kg'], .18)
        self.assertAlmostEqual(result['balances']['l_system']['final_residual_sum_kg'], -1e-4)

    def test_corrupt_signed_integral_is_detected(self):
        rows = fixture(); rows[2]['Mg_in_integral_kg'] += .1
        self.assertFalse(self.analyze(rows)['arithmetic_contract_passes'])

    def test_absolute_transport_is_not_signed_net(self):
        rows = fixture(); rows[2]['Mg_in_abs_integral_kg'] = abs(rows[2]['Mg_in_integral_kg'])
        self.assertFalse(self.analyze(rows)['arithmetic_contract_passes'])

    def test_missing_step_is_rejected(self):
        rows = fixture(); rows[2]['step'] = 3
        with self.assertRaises(ValueError): self.analyze(rows)

    def test_duplicate_step_is_rejected(self):
        rows = fixture(); rows.insert(2, rows[1].copy())
        with self.assertRaises(ValueError): self.analyze(rows)

    def test_wrong_timestep_is_detected(self):
        rows = fixture(); rows[2]['dt_s'] = .3
        self.assertFalse(self.analyze(rows)['arithmetic_contract_passes'])

    def test_nonfinite_accepted_flux_is_rejected(self):
        rows = fixture(); rows[1]['mdot_g_in_kg_s'] = float('nan')
        with self.assertRaises(ValueError): self.analyze(rows)

    def test_initial_unavailable_flux_is_explicit(self):
        rows = fixture(); rows[0]['flux_status'] = 3001; rows[0]['mdot_g_in_kg_s'] = float('nan')
        self.assertTrue(self.analyze(rows)['arithmetic_contract_passes'])

    def test_rejected_flux_status_cannot_enter_ledger(self):
        rows = fixture(); rows[1]['flux_status'] = 3001
        with self.assertRaises(ValueError): self.analyze(rows)

    def test_newton_defect_cannot_be_zeroed(self):
        rows = fixture(); rows[1]['mass_residual_g_pipe_kg'] = 0
        self.assertFalse(self.analyze(rows)['arithmetic_contract_passes'])

    def test_segment_inventory_partition_is_checked(self):
        rows = fixture(); rows[1]['Ml_pipe_kg'] += 1
        self.assertFalse(self.analyze(rows)['arithmetic_contract_passes'])

    def test_negative_inventory_is_not_made_positive(self):
        rows = fixture(); rows[1]['Mg_pipe_kg'] = -.1
        with self.assertRaises(ValueError): self.analyze(rows)

    def test_empty_ledger_is_rejected(self):
        with self.assertRaises(ValueError): self.analyze([])


if __name__ == '__main__':
    unittest.main()
