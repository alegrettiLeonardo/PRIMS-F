import math
import unittest
from apparent_order import apparent_order


class OrderTests(unittest.TestCase):
    def test_unequal_cell_count_ratios_recover_manufactured_order(self):
        h = [57.4/n for n in (51,101,201)]
        for p in (.7,1.,2.4):
            q = [4.2+3*x**p for x in h]
            result = apparent_order(h,q)
            self.assertEqual(result['status'],'POSITIVE_OBSERVED_ORDER')
            self.assertAlmostEqual(result['apparent_order'],p,places=11)

    def test_equal_timestep_halving_recovers_linear_error(self):
        result=apparent_order([.125,.0625,.03125],[10.25,10.125,10.0625])
        self.assertAlmostEqual(result['apparent_order'],1.,places=12)

    def test_error_growth_is_reported_as_nonpositive(self):
        h=[1.,.5,.25];q=[2+x**(-.8) for x in h]
        result=apparent_order(h,q)
        self.assertEqual(result['status'],'NONPOSITIVE_OBSERVED_ORDER')
        self.assertAlmostEqual(result['apparent_order'],-.8,places=12)

    def test_sign_reversal_is_not_promoted_to_convergence(self):
        result=apparent_order([1.,.5,.25],[3.,1.,2.])
        self.assertEqual(result['status'],'NONMONOTONIC_DIFFERENCES')
        self.assertIsNone(result['apparent_order'])
        self.assertAlmostEqual(result['magnitude_order_candidate'],1.,places=12)

    def test_resolution_floor_prevents_order_from_flat_data(self):
        result=apparent_order([1.,.5,.25],[10.0002,10.0001,10.00005],.001)
        self.assertEqual(result['status'],'DIFFERENCE_BELOW_DECLARED_FLOOR')
        self.assertIsNone(result['apparent_order'])

    def test_tiny_sign_reversal_does_not_underflow_away(self):
        result=apparent_order([4.,2.,1.],[3e-200,1e-200,2e-200])
        self.assertEqual(result['status'],'NONMONOTONIC_DIFFERENCES')
        self.assertIsNone(result['apparent_order'])

    def test_invalid_floor_and_overflowing_difference_refused(self):
        for floor in (math.nan,math.inf,-1.):
            with self.assertRaises(ValueError):apparent_order([4.,2.,1.],[3.,2.,1.],floor)
        with self.assertRaises(ValueError):apparent_order([4.,2.,1.],[1e308,-1e308,1e308])

    def test_scale_invariance(self):
        h=[1/51,1/101,1/201];q=[3*x**1.3 for x in h]
        first=apparent_order(h,q)['apparent_order']
        second=apparent_order([100*x for x in h],[-7*x for x in q])['apparent_order']
        self.assertAlmostEqual(first,second,places=12)

    def test_invalid_scales_and_nonfinite_observations_refused(self):
        for h,q in (([1,1,.5],[3,2,1]),([1,.5,0],[3,2,1]),([1,.5,.25],[3,math.nan,1])):
            with self.assertRaises(ValueError):apparent_order(h,q)


if __name__=='__main__':unittest.main()
