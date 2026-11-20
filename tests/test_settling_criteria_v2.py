"""Synthetic qualification of the pre-registered settling/convergence criteria v2."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
P=Path(__file__).resolve().parents[1]/'validation/diagnostics/p0p1c_20261119/settling_criteria.py'
spec=importlib.util.spec_from_file_location('settling_criteria_v2',P);v2=importlib.util.module_from_spec(spec);spec.loader.exec_module(v2)

T=23.4;W=2*np.pi/T
def series(shift=0.0,**mods):
    t=np.arange(0.,320.,.05);x=W*(t-shift)
    m=dict(p_riser_base_Pa=2.66e5+1e4*np.sin(x),p_junction_Pa=2.5e5+9e3*np.sin(x-.3),
           holdup_liquid_base=.6+.2*np.sin(x),holdup_liquid_out=.4+.1*np.sin(x-1),
           mdot_l_base_kg_s=1.28+.5*np.sin(x-.5),mdot_l_out_kg_s=1.28+.8*np.sin(x-1.5),
           mdot_g_base_kg_s=2e-4*np.sin(x-.7),mdot_g_out_kg_s=1e-3+5e-4*np.sin(x-2),
           Ml_kg=71+2*np.sin(x+.4),Ml_pipe_kg=40+1.5*np.sin(x+.4),Ml_riser_kg=31+.5*np.sin(x+.2),
           Mg_kg=.5+.05*np.sin(x-.4),Mg_pipe_kg=.3+.03*np.sin(x-.4),Mg_riser_kg=.2+.02*np.sin(x-.1))
    for k,f in mods.items():m[k]=f(t,m[k])
    edges=[T/4+shift+T*k for k in range(11)]
    cons=dict(time_s=t,Ml_in_integral_kg=1.28*t,Mg_in_integral_kg=1e-3*t)
    return t,m,edges,cons

def run(t,m,edges,cons):
    v2.validate_markers(edges)
    p=m['p_riser_base_Pa'];c0=v2.pressure_c0(np.diff(edges),np.full(10,1e4),
        np.array([v2.window_mean(t,p,l,r) for l,r in zip(edges[:-1],edges[1:])]))
    return v2.settle(t,m,edges,v2.feed_rates(cons,edges),c0)

class Settling(unittest.TestCase):
    def test_periodic_limit_cycle_settles(self):
        r=run(*series());self.assertTrue(r['settled'],r['failed_criteria']);self.assertEqual(r['missing_channels'],[])
    def test_drifting_inventory_fails_despite_settled_pressure(self):
        t,m,e,c=series(Ml_riser_kg=lambda t,y:y+.003*t)  # 0.7 kg over the window, ~70% of its range
        r=run(t,m,e,c)
        self.assertTrue(r['pressure_C0']['passed'])
        self.assertFalse(r['settled']);self.assertIn('Ml_riser_kg:drift',r['failed_criteria'])
    def test_slow_drift_small_relative_to_mean_is_still_caught(self):
        # 0.3% of the 31 kg mean: invisible to a mean-normalized 0.5% rule.
        t,m,e,c=series(Ml_riser_kg=lambda t,y:y+.093*(t-t[0])/(t[-1]-t[0]))
        r=run(t,m,e,c);self.assertIn('Ml_riser_kg:drift',r['failed_criteria'])
    def test_varying_waveform_with_equal_mean_and_range_fails_shape(self):
        def alt(t,y):
            k=np.floor((t-T/4)/T);s=np.where(k%2==0,1.,-1.)
            return y+.04*s*np.sin(2*W*t)
        r=run(*series(holdup_liquid_base=alt))
        self.assertTrue(r['pressure_C0']['passed'])
        self.assertIn('holdup_liquid_base:shape',r['failed_criteria'])
        self.assertNotIn('holdup_liquid_base:drift',r['failed_criteria'])
    def test_near_zero_mean_channel_uses_range_not_mean(self):
        r=run(*series());g=r['channels']['mdot_g_base_kg_s']
        self.assertLess(abs(g['mean']),1e-6);self.assertEqual(g['scale'],g['mean_range'])
    def test_constant_channel_uses_dimensioned_floor(self):
        r=run(*series(mdot_g_base_kg_s=lambda t,y:0*y));g=r['channels']['mdot_g_base_kg_s']
        self.assertEqual(g['mean_range'],0.);self.assertAlmostEqual(g['scale'],1e-4*1e-3,places=15);self.assertTrue(g['passed'])
    def test_zero_feed_zero_range_is_not_evaluable(self):
        t,m,e,c=series(mdot_g_base_kg_s=lambda t,y:0*y);c['Mg_in_integral_kg']=0*t
        r=run(t,m,e,c);self.assertIn('mdot_g_base_kg_s',r['missing_channels']);self.assertFalse(r['settled'])
    def test_missing_channel_is_not_settled(self):
        t,m,e,c=series();del m['Mg_riser_kg'];r=run(t,m,e,c)
        self.assertFalse(r['settled']);self.assertIn('Mg_riser_kg',r['missing_channels'])
    def test_eleven_markers_required(self):
        t,m,e,c=series()
        with self.assertRaises(ValueError):run(t,m,e[:10],c)
        with self.assertRaises(ValueError):run(t,m,e[:5]+[e[4]]+e[6:],c)

class Convergence(unittest.TestCase):
    def test_identical_runs_converge(self):
        a=run(*series());r=v2.compare(a,a);self.assertTrue(r['converged_v2_channels'])
    def test_phase_alignment_cannot_hide_phase_error(self):
        a=run(*series())
        # Holdup shifted by 5% of the period relative to the pressure markers.
        b=run(*series(holdup_liquid_base=lambda t,y:.6+.2*np.sin(W*(t-.05*T))))
        r=v2.compare(a,b);ch=r['channels']['holdup_liquid_base']
        self.assertFalse(ch['passes']['shape_phase_fixed'])
        self.assertAlmostEqual(abs(ch['phase_error_fraction_diagnostic']),.05,delta=1/256)
        self.assertFalse(r['converged_v2_channels'])
    def test_inventory_level_difference_fails_even_with_equal_pressure(self):
        a=run(*series());b=run(*series(Ml_kg=lambda t,y:y+1.0))
        r=v2.compare(a,b);self.assertIn('Ml_kg:level',r['failed_criteria'])
    def test_unsettled_runs_never_converge(self):
        a=run(*series());b=run(*series(Ml_riser_kg=lambda t,y:y+.003*t))
        r=v2.compare(a,b);self.assertFalse(r['both_settled']);self.assertFalse(r['converged_v2_channels'])

if __name__=='__main__':unittest.main()
