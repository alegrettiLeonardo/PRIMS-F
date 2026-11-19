import json
from offline_probe import Probe,InvalidState,ROOT
p=Probe('D_SS1',3)
result={'mode':'Rusanov only at MATCHED, original Roe spectral radius, no multiplier', 'states':{}}
for name,x in [('accepted_old',p.old),('linearization',p.trial)]:
    d=p.interface(x)
    item={'gas_mass_flux_kg_s':float(d[1]*d[45]*d[46]*d[47]),'jg_m_s':float(p.margins(x)[0]),
          'minimum_phase_masses_normalized':x.reshape(-1,3)[:,:2].min(axis=0).tolist()}
    try:
        item['residual_norm']=p.evaluate(x,False)[2]
        item['status']=0
    except InvalidState as error:
        item.update(error.details)
    result['states'][name]=item
(ROOT/'rusanov_interface.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
