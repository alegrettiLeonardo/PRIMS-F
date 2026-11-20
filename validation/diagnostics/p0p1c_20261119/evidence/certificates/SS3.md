# Certificado de bloqueio — SS3

Execução: `/home/leonardo/Documentos/PRIMS-F-p0p1c-results-20261119/ss3/SS3_CENTER`

* Bloqueio em: **NEWTON_PROPOSAL**; classes: **NO_ROOT_IN_DECLARED_BRANCH, CONSTITUTIVE_DOMAIN_LIMIT**
* t = 97.1318521911 s; último dt aceito 2.214e-07 s; menor tentado 1.329e-07 s
* Local: {'index': 69, 's_m': 61.2131, 'inclination_deg': 51.6372, 'width_m': 0.246005}
* Fronteira do raio: {'lambda_admissible': 3.6817117703119605e-06, 'lambda_refused': 3.681711770311961e-06, 'status': 1017, 'index': 69, 'stage': 'CELL', 'norm_ratio_at_edge': 0.9999963, 'decoded': {'site': 'cell primitive recovery', 'primitive_status': 7, 'primitive_name': 'NO_PHYSICAL_ROOT'}}
* Jacobiano: {'lapack_info': 0, 'norm1': 1.358046, 'condition_1norm_estimate': 1.844256, 'fd_max_relative_error_local_columns': 1.877685e-09, 'fd_columns_without_admissible_stencil': 0}
* Raio: {'admissible_samples': 210, 'best_norm': 2.242448e-08, 'best_lambda': 3.667491e-06, 'best_norm_ratio': 0.9999963}
* Margens: {'status': 0, 'junction_raw_gas_superficial_nd': 0.2269349, 'min_riser_drift_margin_m_s': 9.574286e-13, 'at_cell': 69}

## Evidências
* Jacobian assembled at the trial, LU nonsingular, local FD agreement 1.9e-09: not a Jacobian failure
* admissible part of the Newton ray reduces the residual only to 0.999996: no admissible descent, not a line-search weakness
* primitive recovery status 7: analytic exclusion over the whole co-current interval
* root of the trial lies at the co-current endpoint (u_l=0); the proposal requires u_l<0

_absence of roots along local Newton proposals is not a proof that no global solution exists; the model, not the physics, is certified here_
