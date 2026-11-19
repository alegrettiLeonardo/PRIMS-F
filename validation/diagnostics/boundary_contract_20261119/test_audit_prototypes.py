"""Regression tests for isolated prototypes; production source files are immutable."""
from pathlib import Path
import argparse,unittest,sys
import numpy as np
from audit_local import Audit,fd
BASE=Path('/home/leonardoab/Documentos/PRIMS-F-boundary-results-20261119')
class Prototypes(unittest.TestCase):
 def audit(self,v,debug=False):return Audit(BASE/('variants_v2_debug' if debug else 'variants_v2')/v/'audit.so')
 def test_volume_contract_and_physical_flux(self):
  for debug in [False,True]:
   a=self.audit('volume',debug)
   for p in [1.,2.1,3.]:
    u=a.state(p=p);b=a.boundary(u);self.assertEqual(b['status'],0)
    rho=1+(b['p']-1)/a.par[12]**2
    self.assertLess(abs(b['f'][0]/rho/a.par[19]-1),1e-10)
    point=a.point(b['bc']);self.assertEqual(point['status'],0)
    self.assertLess(np.max(abs(point['f']-b['f'])/np.maximum(1,abs(b['f']))),1e-10)
 def test_partition_invariance(self):
  for debug in [False,True]:
   a=self.audit('partition',debug)
   u=np.array([a.state(p=2.1+.001*i,alpha=.25+.01*i,mom=.6+.003*i) for i in range(6)])
   first=a.assemble(u,npip=2);last=a.assemble(u,npip=4)
   self.assertEqual([first['status'],last['status']],[0,0]);self.assertLess(np.max(abs(first['r']-last['r'])),1e-12)
 def test_full_jacobian_with_independent_stencil(self):
  for debug in [False,True]:
   for v in ['standard','volume','partition']:
    a=self.audit(v,debug);u=a.state();old=u[None,:];r=a.assemble(old,npip=1,jac=1);self.assertEqual(r['status'],0)
    for rel in [1e-3,1e-4,1e-5]:
     independent=fd(lambda x:a.assemble(x.reshape(1,3),old=old,npip=1)['r'].ravel(order='F'),u,rel)
     self.assertLess(np.max(abs(independent-r['j'])/np.maximum(1,abs(independent))),1e-5)
 def test_zero_flow_normalization_is_undefined(self):
  from analyze_audit import normalized
  self.assertIsNone(normalized(0.,0.));self.assertIsNone(normalized(1.,0.));self.assertEqual(normalized(2.,4.),.5)
 def test_absent_phase_has_explicit_derivative_status(self):
  a=self.audit('standard');a.par[18]=0;a.setup();u=a.state(alpha=0,mom=.5)[None,:]
  flux=a.boundary(u[0]);self.assertEqual(flux['status'],0);self.assertEqual(flux['bc'][1],0.)
  op=a.assemble(u,npip=1,jac=1);self.assertEqual(op['status'],1031);self.assertEqual(op['index'],1)
 def test_uniform_stationary_state(self):
  a=self.audit('standard');a.par[8]=0;a.par[17:20]=0;a.setup()
  u=np.repeat(a.state(alpha=.25,mom=0)[None,:],4,axis=0);op=a.assemble(u,npip=2,jac=1)
  self.assertEqual(op['status'],0);self.assertLess(np.max(abs(op['r'])),1e-12)
if __name__=='__main__':unittest.main()
