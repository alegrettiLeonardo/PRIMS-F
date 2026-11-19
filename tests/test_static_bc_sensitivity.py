from decimal import Decimal, localcontext
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_static_bc_sensitivity import offsets, validate_nominal


def nominal():
    return dict(nodes=51.,n_riser=51.,dt_max_s=.125,time_dim=3600.,tol_res=1e-8,
                tol_inc=1e-8,ps=206000.,jg=1.16,jl=.547,tgas=293.15,
                gas_flow_basis='STANDARD',p_standard=101300.,t_standard=293.15,
                rg=287.,rho_l0=997.05,mul=545.38635,mug=1.396675899947762,
                pipe_regime_mode='FIXED_STRATIFIED',closure_riser=3.,
                gas_eos='BAROTROPIC',dissipation_mode='LAGGED',
                junction_model='MATCHED',flux_scheme='ROE')


class StaticBoundarySensitivityTest(unittest.TestCase):
    def test_six_static_axis_pairs_and_no_forcing(self):
        values=nominal();saved=dict(values);cases=offsets(values)
        self.assertEqual(values,saved)
        self.assertEqual(len(cases),6)
        self.assertEqual([axis for _,axis,_ in cases],['ps','ps','jg','jg','jl','jl'])
        for _,axis,changes in cases:
            self.assertEqual(set(changes),{'ps'} if axis=='ps' else {axis,'mug' if axis=='jg' else 'mul'})
        self.assertEqual([c['ps'] for _,a,c in cases if a=='ps'],[201000.,211000.])

    def test_si_overrides_against_decimal_reference(self):
        with localcontext() as context:
            context.prec=50
            for _,axis,changes in offsets(nominal()):
                if axis=='jg':
                    reference=Decimal('101300')*Decimal(str(changes[axis]))/(Decimal('287')*Decimal('293.15'))
                    self.assertAlmostEqual(changes['mug'],float(reference),places=14)
                elif axis=='jl':
                    reference=Decimal('997.05')*Decimal(str(changes[axis]))
                    self.assertAlmostEqual(changes['mul'],float(reference),places=12)

    def test_wrong_flow_basis_or_nominal_is_rejected(self):
        for key,value in [('gas_flow_basis','ACTUAL'),('ps',201000.),('time_dim',2400.),('tol_res',1e-10)]:
            values=nominal();values[key]=value
            with self.assertRaises(ValueError):validate_nominal(values)

    def test_inconsistent_mass_override_rejected(self):
        for key in ('mul','mug'):
            values=nominal();values[key]*=1.01
            with self.assertRaises(ValueError):offsets(values)


if __name__=='__main__':unittest.main()
