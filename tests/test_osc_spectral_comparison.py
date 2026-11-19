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


if __name__=="__main__":
    unittest.main()
