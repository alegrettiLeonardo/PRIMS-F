import math,unittest
from analyze_fixed_trace_holdup import calculate

class FixedTraceTests(unittest.TestCase):
    def fixture(self):
        times=[i/1000 for i in range(11001)]
        return {'time_s':times,'alpha_riser_trace':[.4-.2*math.sin(2*math.pi*t) for t in times]}

    def test_known_periodic_holdup_without_amplitude_rescaling(self):
        out=calculate(self.fixture(),list(range(1,12)))
        self.assertTrue(out['auxiliary_accommodated'])
        self.assertAlmostEqual(out['cycle_mean_statistics']['mean'],.6,places=12)
        self.assertAlmostEqual(out['cycle_mean_statistics']['peak_to_peak'],.4,places=12)
        self.assertAlmostEqual(max(out['waveform']['mean']),.8,places=12)
        self.assertAlmostEqual(min(out['waveform']['mean']),.4,places=12)

    def test_missing_and_invalid_fraction_refused_inside_window(self):
        for value in (None,math.nan,1.01,-.01):
            data=self.fixture();data['alpha_riser_trace'][2000]=value
            with self.assertRaises(ValueError):calculate(data,list(range(1,12)))

    def test_missing_outside_window_does_not_pollute_selected_cycles(self):
        data=self.fixture();data['alpha_riser_trace'][0]=None
        self.assertTrue(calculate(data,list(range(1,12)))['auxiliary_accommodated'])

    def test_extrapolation_and_nonmonotonic_times_refused(self):
        data=self.fixture()
        with self.assertRaises(ValueError):calculate(data,list(range(2,13)))
        data['time_s'][2000]=data['time_s'][1999]
        with self.assertRaises(ValueError):calculate(data,list(range(1,12)))

if __name__=='__main__':unittest.main()
