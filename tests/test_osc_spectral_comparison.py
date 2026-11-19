import unittest
import numpy as np
from summarize_osc_spectra import spectral_peaks


class SpectralComparisonTests(unittest.TestCase):
    def test_resolves_a_fundamental_and_harmonic_with_finite_bin_width(self):
        time=np.arange(0.,300.,.125)
        signal=4+np.sin(2*np.pi*time*.04)+.4*np.sin(2*np.pi*time*.08)
        result=spectral_peaks(time,signal)
        peaks=result["three_largest_local_peaks"]
        self.assertAlmostEqual(peaks[0]["frequency_bin_Hz"],.04,places=10)
        self.assertAlmostEqual(peaks[1]["frequency_bin_Hz"],.08,places=10)
        self.assertAlmostEqual(result["ratio_two_strongest_frequency_bins"],2.,places=10)
        self.assertAlmostEqual(result["frequency_resolution_Hz"],1/300,places=10)
        self.assertAlmostEqual(peaks[1]["relative_power"],.16,places=4)


    def test_constant_pressure_has_no_invented_spectral_peak(self):
        time=np.arange(0.,300.,.125)
        result=spectral_peaks(time,np.full(len(time),206000.))
        self.assertEqual(result["three_largest_local_peaks"],[])
        self.assertIsNone(result["ratio_two_strongest_frequency_bins"])
        self.assertFalse(result["window_excursion_at_least_100Pa"])

    def test_normalized_small_signal_is_not_qualified_as_a_pressure_cycle(self):
        time=np.arange(0.,300.,.125)
        signal=206000.+.2*np.sin(2*np.pi*time*.04)
        result=spectral_peaks(time,signal)
        self.assertFalse(result["window_excursion_at_least_100Pa"])
        self.assertAlmostEqual(result["window_excursion_Pa"],.4)
        self.assertAlmostEqual(result["three_largest_local_peaks"][0]["frequency_bin_Hz"],.04)
        self.assertIn("do not qualify",result["peaks_interpretation"])


if __name__=="__main__":
    unittest.main()
