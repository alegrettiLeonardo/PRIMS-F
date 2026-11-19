"""Regression for the observed one-step loss in continuation postprocessing."""
import importlib.util,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'validation/diagnostics/p0p1_20261119/step_metadata.py'
spec=importlib.util.spec_from_file_location('step_metadata_test',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class StepMetadata(unittest.TestCase):
 def test_actual_steps_file_has_no_initial_row(self):
  rows=[dict(step=1,dt_s=.001),dict(step=2,dt_s=.002)]
  s=m.statistics(rows)
  self.assertEqual(s['accepted'],2);self.assertEqual(s['min'],.001)
  self.assertAlmostEqual(s['sum_duration_s'],.003)
 def test_initial_ledger_row_is_not_a_step(self):
  rows=[dict(step=0,dt_s=0),dict(step=1,dt_s=.001)]
  self.assertEqual(m.statistics(rows)['accepted'],1)
 def test_first_step_of_each_continuation_is_retained(self):
  first=[dict(step=1,dt_s=.001),dict(step=2,dt_s=.002)]
  second=[dict(step=1,dt_s=.125)]
  combined=m.accepted_rows(first)+m.accepted_rows(second)
  self.assertEqual(len(combined),3)
  self.assertAlmostEqual(m.statistics(combined)['sum_duration_s'],.128)
