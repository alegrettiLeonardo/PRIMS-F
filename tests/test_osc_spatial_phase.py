import unittest
import numpy as np
from compare_osc_spatial_phase import crossings, sample_positions, phase_average, compare_fields


class SpatialPhaseTest(unittest.TestCase):
    def test_affine_physical_profile_reproduced_on_unequal_grids(self):
        target = np.linspace(.2, .8, 7)
        for n in (11, 32, 101):
            x = np.linspace(.1, .9, n)
            value = np.stack((3*x-2, -4*x+9), axis=1)
            np.testing.assert_allclose(sample_positions(x, value, target),
                                       np.stack((3*target-2, -4*target+9), axis=1), atol=2e-15)
        with self.assertRaises(ValueError):
            sample_positions(x, value, [0, .5])

    def test_translating_isovoid_is_geometric_level_not_material_identity(self):
        x = np.array([0., .2, .5, .8, 1.])
        for t in (0, .1, .2):
            result = crossings(x, x-t, .5)
            self.assertEqual(len(result), 1)
            self.assertAlmostEqual(result[0][0], .5+t)
            self.assertEqual(result[0][1], 'increasing')
        self.assertEqual(crossings(x, 1-x, .5), [(.5, 'decreasing')])
        self.assertEqual(crossings(x, np.ones(5), .5), [])

    def test_phase_lag_and_spatial_profile_are_preserved(self):
        time = np.linspace(0, 10, 10001)
        position = np.array([.2, .5, .8])
        signal = np.sin(2*np.pi*(time[:, None]-.25))+position[None, :]
        fields = np.repeat(signal[:, :, None], 8, axis=2)
        mean, scatter = phase_average(time, fields, list(range(11)))
        self.assertLess(scatter.max(), 1e-14)
        self.assertEqual(np.argmax(mean[:, 0, 0]), 128)
        a = {'pipe': mean, 'riser': mean}
        b = {'pipe': mean+2, 'riser': mean+2}
        result = compare_fields(a, b)
        self.assertAlmostEqual(result['pipe']['p_Pa']['rms_absolute'], 2.)
        self.assertAlmostEqual(result['riser']['jg_m_s']['mean_difference'], 2.)


if __name__ == '__main__':
    unittest.main()
