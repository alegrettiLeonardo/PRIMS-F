import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from assess_burr15_primary_bounds import assess, geometry, POINTS


class PrimaryBoundsTest(unittest.TestCase):
    def test_independent_circle_half_and_third(self):
        d = 0.06
        for phi in (math.pi/2, math.pi/3, 2*math.pi/3):
            alpha = 1-(phi-math.sin(phi)*math.cos(phi))/math.pi
            h, sl = geometry(alpha, d)
            self.assertAlmostEqual(h, d*(1-math.cos(phi))/2, places=15)
            self.assertAlmostEqual(sl, d*phi, places=15)

    def test_contradiction_survives_reported_uncertainty_envelope(self):
        for point in POINTS:
            out = assess(point)
            self.assertTrue(out['literal_detector_loss_for_entire_envelope'])
            self.assertGreater(out['signed_margin_lower_m_s'], 0.34)
            self.assertLess(out['jg_envelope_m_s'][1], 1.0)

    def test_upper_bound_contains_dense_interior_states(self):
        for point in POINTS:
            out = assess(point)
            low, high = out['alpha_envelope']
            for k in range(101):
                a = low + (high-low)*k/100
                h, sl = geometry(a, 0.06)
                critical = (1-h/0.06)*math.sqrt(9.81*math.pi*0.06**2/4*a/sl)
                self.assertLessEqual(critical, out['burr15_critical_upper_m_s']*(1+1e-14))
                self.assertGreaterEqual(out['jg_envelope_m_s'][0]/a,
                                        out['ug_lower_m_s']*(1-1e-14))

    def test_absent_phases_not_extrapolated(self):
        for alpha in (-0.1, 0, 1, 1.1, math.nan):
            with self.assertRaises(ValueError):
                geometry(alpha, 0.06)


if __name__ == '__main__':
    unittest.main()
