# Certificado de bloqueio — OSC_ISOTHERMAL

Execução: `/home/leonardo/Documentos/PRIMS-F-p0p1c-results-20261119/iso/B_ISO_MASS`

* Bloqueio em: **NEWTON_PROPOSAL**; classes: **CONSTITUTIVE_DOMAIN_LIMIT, UNSUPPORTED_PHYSICAL_BOUNDARY**
* t = 223.6445229171 s; último dt aceito 1.653e-07 s; menor tentado 1.984e-07 s
* Local: {'index': 53, 's_m': 56.8373, 'inclination_deg': -2.0, 'width_m': 1.12549}
* Fronteira do raio: {'lambda_admissible': 3.3698803140518686e-08, 'lambda_refused': 3.369880314051869e-08, 'status': 3012, 'index': 53, 'stage': 'FACE', 'norm_ratio_at_edge': 1.0, 'decoded': {'site': 'junction riser closure', 'closure_status': 12}}
* Jacobiano: {'lapack_info': 0, 'norm1': 8.032173, 'condition_1norm_estimate': 63.93681, 'fd_max_relative_error_local_columns': 2.753465e-06, 'fd_columns_without_admissible_stencil': 0}
* Raio: {'admissible_samples': 212, 'best_norm': 2.901624e-07, 'best_lambda': 0.0002262828, 'best_norm_ratio': 0.9997737}
* Margens: {'status': 0, 'junction_raw_gas_superficial_nd': 1.80504e-13, 'min_riser_drift_margin_m_s': 17.05364, 'at_cell': 103}

## Evidências
* Jacobian assembled at the trial, LU nonsingular, local FD agreement 2.8e-06: not a Jacobian failure
* admissible part of the Newton ray reduces the residual only to 0.999774: no admissible descent, not a line-search weakness
* first refused point: CLOSURE_DOMAIN of the junction riser closure
* junction gas flux 5.81e-13 -> -9.27e-13 kg/m2/s with face riser holdup 1.99e-13: required mixed->pure liquid transition of the face trace while both cells stay two-phase

_absence of roots along local Newton proposals is not a proof that no global solution exists; the model, not the physics, is certified here_
