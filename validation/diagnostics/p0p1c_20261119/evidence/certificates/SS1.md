# Certificado de bloqueio — SS1

Execução: `/home/leonardo/Documentos/PRIMS-F-p0p1c-results-20261119/jt2/SS1_CENTER`

* Bloqueio em: **NEWTON_PROPOSAL**; classes: **CONSTITUTIVE_DOMAIN_LIMIT, UNSUPPORTED_PHYSICAL_BOUNDARY**
* t = 78.5766182805 s; último dt aceito 1.538e-07 s; menor tentado 1.845e-07 s
* Local: {'index': 53, 's_m': 56.8373, 'inclination_deg': -2.0, 'width_m': 1.12549}
* Fronteira do raio: {'lambda_admissible': 1.5880661467349782e-09, 'lambda_refused': 1.5880661467349784e-09, 'status': 3012, 'index': 53, 'stage': 'FACE', 'norm_ratio_at_edge': 1.0, 'decoded': {'site': 'junction riser closure', 'closure_status': 12}}
* Jacobiano: {'lapack_info': 0, 'norm1': 7.209976, 'condition_1norm_estimate': 41.17084, 'fd_max_relative_error_local_columns': 8.834594e-07, 'fd_columns_without_admissible_stencil': 0}
* Raio: {'admissible_samples': 161, 'best_norm': 1.420781e-07, 'best_lambda': 7.647589e-05, 'best_norm_ratio': 0.9999235}
* Margens: {'status': 0, 'junction_raw_gas_superficial_nd': 2.859629e-15, 'min_riser_drift_margin_m_s': 370.4085, 'at_cell': 85}

## Evidências
* Jacobian assembled at the trial, LU nonsingular, local FD agreement 8.8e-07: not a Jacobian failure
* admissible part of the Newton ray reduces the residual only to 0.999923: no admissible descent, not a line-search weakness
* first refused point: CLOSURE_DOMAIN of the junction riser closure
* junction gas flux 6.85e-15 -> -2.32e-14 kg/m2/s with face riser holdup 3.20e-15: required mixed->pure liquid transition of the face trace while both cells stay two-phase
* both traces still carry gas upward (7.115e-02, 4.038e-02 kg/m2/s): the vanishing face flux comes from the Roe dissipation (diss_Gg=-5.576e-02) at the junction holdup jump; a PHYSICAL blockage is NOT established by these states

_absence of roots along local Newton proposals is not a proof that no global solution exists; the model, not the physics, is certified here_
