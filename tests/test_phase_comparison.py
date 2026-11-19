import unittest
import numpy as np
from compare_osc_phase import phase_waveform,compare


class PhaseComparisonTests(unittest.TestCase):
    def test_known_phase_and_amplitude_recovered_without_optimized_shift(self):
        time=np.linspace(0,100,10001)
        signal=20+3*np.cos(2*np.pi*(time/10-.25))
        wave=phase_waveform(time,signal,list(range(0,101,10)))
        self.assertAlmostEqual(wave["excursion_of_phase_mean"],6,places=5)
        self.assertAlmostEqual(wave["phase_at_maximum_relative_to_pressure_peak"],.25)
        self.assertLess(wave["phase_cycle_scatter_rms"],1e-12)

    def test_flat_channel_has_no_invented_phase_marker(self):
        wave=phase_waveform(np.arange(101),np.ones(101),list(range(0,101,10)))
        self.assertEqual(wave["excursion_of_phase_mean"],0)
        self.assertIsNone(wave["phase_at_maximum_relative_to_pressure_peak"])

    def test_no_cycle_comparison_reports_status_without_waveforms(self):
        left={"run":"flat","freeze":"a","binary_sha256":"b","controls":{},"status":"NO_CYCLES",
              "cycles":{"available":False}}
        result=compare(left,dict(left,run="flat2"))
        self.assertFalse(result["both_settled"])
        self.assertNotIn("mean_period_s",result)


if __name__=="__main__":
    unittest.main()
