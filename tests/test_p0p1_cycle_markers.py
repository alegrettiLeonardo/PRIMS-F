"""Regress the eleven-marker API consumed by the new ten-cycle qualification."""
import importlib.util
from pathlib import Path
import unittest
P=Path(__file__).resolve().parents[1]/'validation/diagnostics/p0p1_20261119/analyze_refinement.py'
spec=importlib.util.spec_from_file_location('refinement_markers',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Markers(unittest.TestCase):
 def test_terminal_marker_is_not_an_extra_cycle(self):
  summary=dict(peak_records=[dict(peak_time_s=20.*i)for i in range(11)],window_end_s=200.)
  edges=m.cycle_markers(summary)
  self.assertEqual(len(edges),11)
  self.assertEqual([r-l for l,r in zip(edges[:-1],edges[1:])],[20.]*10)
 def test_ten_markers_are_only_nine_complete_cycles(self):
  with self.assertRaises(ValueError):m.cycle_markers(dict(peak_records=[dict(peak_time_s=20.*i)for i in range(10)],window_end_s=180.))
