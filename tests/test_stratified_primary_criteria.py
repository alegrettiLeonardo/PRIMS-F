import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_stratified_primary_criteria import compare, td1977_factor
from assess_burr15_primary_bounds import POINTS, assess


class PrimaryCriteriaTest(unittest.TestCase):
    def test_half_full_analytic_geometry_and_diameter_scaling(self):
        k=td1977_factor(.5,.06,9.81)
        self.assertAlmostEqual(k,9.81*math.pi*.06/32,places=14)
        self.assertAlmostEqual(td1977_factor(.5,.24,9.81),4*k,places=14)

    def test_conditional_threshold_against_dense_rectangle(self):
        for point in POINTS:
            report=compare(point);bounds=assess(point)
            ratio=report['td1977_no_loss_entire_envelope_if_rhol_over_rhog_strictly_greater_than']
            a,b=bounds['alpha_envelope']
            for i in range(101):
                alpha=a+(b-a)*i/100
                ug=bounds['jg_envelope_m_s'][1]/alpha
                self.assertLessEqual(ug**2,(ratio-1)*td1977_factor(alpha,.06,9.81)*(1+1e-14))

    def test_no_unreferenced_angle_extension(self):
        for angle in (-.01,.01,math.nan):
            with self.assertRaises(ValueError):td1977_factor(.7,.06,9.81,angle)

if __name__=='__main__':unittest.main()
