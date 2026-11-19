import unittest
from analyze_cycle_conservation import cycle_balances


def manufactured():
    result=[]
    for t in (0.,.7,1.9,3.1):
        row={'time_s':t}
        for p in ('l','g'):
            flows={'in':1.,'out':.5,'base_left':.8,'base_right':.8} if p=='l' else {'in':-.1,'out':-.3,'base_left':-.2,'base_right':-.2}
            pipe=10.+(flows['in']-flows['base_left']+1e-4)*t
            riser=20.+(flows['base_right']-flows['out']-2e-4)*t
            for r,value,slope in (('',pipe+riser,-1e-4),('_pipe',pipe,1e-4),('_riser',riser,-2e-4)):
                row[f'M{p}{r}_kg']=value;row[f'mass_residual_sum_{p}{r}_kg']=slope*t
            for face,flow in flows.items():
                row[f'M{p}_{face}_integral_kg']=flow*t;row[f'M{p}_{face}_abs_integral_kg']=abs(flow)*t
        result.append(row)
    return result


class CycleBalanceTests(unittest.TestCase):
    def test_nonuniform_steps_and_partial_step_boundaries(self):
        cycles=cycle_balances(iter(manufactured()),[.1,1.,2.4,3.])
        self.assertEqual(len(cycles),3)
        for cycle in cycles:
            dt=cycle['period_s']
            for p in ('l','g'):
                for r,slope in (('',-1e-4),('_pipe',1e-4),('_riser',-2e-4)):
                    metric=cycle['balances'][p+(r or '_system')]
                    self.assertAlmostEqual(metric['final_balance_kg'],slope*dt,places=12)
                    self.assertAlmostEqual(metric['raw_balance_vs_residual_difference_kg'],0.,places=12)
                    self.assertFalse(metric['passes_inventory_1e6'])

    def test_countercurrent_absolute_transport_retains_signs(self):
        cycle=cycle_balances(iter(manufactured()),[0.,3.])[0]
        self.assertAlmostEqual(cycle['balances']['g_system']['absolute_processed_mass_kg'],1.2,places=12)
        self.assertAlmostEqual(cycle['balances']['g_pipe']['absolute_processed_mass_kg'],.9,places=12)

    def test_missing_coverage_and_reversed_edges_rejected(self):
        for edges in ((-.1,1.),(0.,4.),(1.,.5)):
            with self.assertRaises(ValueError):cycle_balances(iter(manufactured()),edges)

    def test_duplicate_accepted_time_rejected(self):
        rows=manufactured();rows.insert(2,rows[1].copy())
        with self.assertRaises(ValueError):cycle_balances(iter(rows),[0.,3.])


if __name__=='__main__':unittest.main()
