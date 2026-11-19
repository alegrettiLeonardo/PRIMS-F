"""Compare interface mass properties, without inventing undefined momentum traces."""
import json
from offline_probe import Probe, InvalidState, ROOT, np

p=Probe('D_SS1',1)
width=np.loadtxt(ROOT/'D_SS1/last_rejected_state.dat',skiprows=2)[:,2]
wu,rho,pref=p.meta[3:]
pip=pref/(rho*wu**2)
cg=343*np.sqrt(rho/pref)
dt_nd=p.meta[1]/(12.546253993197080/wu)
report={'case':'D_SS1','source_operator':'69813ee','states':{},
        'hll_source':'Harten, Lax, van Leer1983, p50, Eqs3.13–3.15',
        'wave_bound_limit':'Endpoint eigenvalue/phase-speed envelope; no proof that it bounds every nonlinear wave.',
        'positivity_limit':'Mapped HLL intermediate and explicit mass update are checked separately from trace domain.'}
for label,x in [('accepted_old',p.old),('linearization',p.trial)]:
    data=p.interface(x)
    raw=data[:3]
    hll=data[36:39]
    ps=data[3]
    kin=raw[2]-pip*ps
    hllps=(hll[2]-kin)/pip
    item={'roe_gas_mass_flux_kg_s':float(raw[1]*rho*wu*data[47]),
          'roe_jg_m_s':float(raw[1]/(ps/cg**2)*wu),
          'hll_mass_fluxes_kg_s':(hll[:2]*rho*wu*data[47]).tolist(),
          'hll_pressure_reconstruction_Pa':float(hllps*pref),
          'hll_jg_m_s':float(hll[1]/(hllps/cg**2)*wu),
          'hll_mapped_intermediate_phase_mass_density_kg_m3':(data[39:41]*rho).tolist(),
          'mapped_vs_original_right_phase_mass_ratio':(data[25:27]/data[42:44]).tolist(),
          'wave_speeds_m_s':(data[34:36]*wu).tolist(),
          'hll_trace_admissible_bhagwat':bool(hll[1]>=0),
          'physical_cell_mass_min_kg_m3':(x.reshape(-1,3)[:,:2].min(axis=0)*rho).tolist()}
    if label=='accepted_old':
        residual=p.evaluate(x,False)[0].reshape(-1,3)
        alternate=residual.copy()
        delta=dt_nd*(hll[:2]-raw[:2])
        alternate[50,:2]+=delta/width[50]
        alternate[51,:2]-=delta/width[51]
        updated=x.reshape(-1,3)-alternate
        item['explicit_mass_update_min_kg_m3']=(updated[:,:2].min(axis=0)*rho).tolist()
        item['interface_mass_cancellation_nd']=np.sum((alternate-residual)[:,:2]*width[:,None],axis=0).tolist()
        item['momentum_update_defined']=False
    report['states'][label]=item
try:
    # This is the existing global RUSANOV option, with no spectral multiplier.
    p=Probe('D_SS1',2)
    report['existing_rusanov_prepare_status']=0
except InvalidState as error:
    report['existing_rusanov_prepare_failure']=error.details
(ROOT/'flux_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
