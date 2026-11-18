# Final frozen closure comparison

Source `490936a676c2a86f1348f0c867bda65685c27823` includes the derivative-tail correction and stricter CCFL singularity/error propagation. Both closures, the stable control, all unstable attempts and the independent mesh/time refinements were rerun on this single source/binary version.

| Case | Closure | Cells/section | dt max [s] | End [s] | Result | Solver status |
|---|---:|---:|---:|---:|---|---:|
| D_OSC | 3 | 51 | 0.125 | 2400.000000000 | COMPLETED | 0 |
| D_OSC | 5 | 51 | 0.125 | 93.100004084 | FAILED | 6000 |
| D_SS1 | 3 | 51 | 0.125 | 78.571788477 | FAILED | 6000 |
| D_SS1 | 5 | 51 | 0.125 | 52.795265494 | FAILED | 3061 |
| D_SS3 | 3 | 51 | 0.125 | 97.131852339 | FAILED | 6000 |
| D_SS3 | 5 | 51 | 0.125 | 42.182009549 | FAILED | 6000 |
| D_Stable | 3 | 51 | 0.125 | 2400.000000000 | COMPLETED | 0 |
| D_Stable | 5 | 51 | 0.125 | 2400.000000000 | COMPLETED | 0 |
| D_OSC | 3 | 51 | 0.0625 | 2400.000000000 | COMPLETED | 0 |
| D_OSC | 5 | 51 | 0.0625 | 77.364993950 | FAILED | 6000 |
| D_OSC | 3 | 101 | 0.125 | 2400.000000000 | COMPLETED | 0 |
| D_OSC | 5 | 101 | 0.125 | 94.245285780 | FAILED | 6000 |
| D_OSC | 3 | 51 | 0.125 | 3600.000000000 | COMPLETED | 0 |

Closure 3: Bhagwat. Closure 5: Chexal1992 AUTO. Requested horizons and all raw failures remain explicit. A truncated failed run is not assigned a full-window period.

The SS1 MATCHED mapping still proposes a negative horizontal gas trace while its last accepted/trial cell states remain admissible. Chexal returns DOMAIN61 for this trace; no projection was introduced. In Chexal SS3/OSC, the last saved states retain tiny positive liquid velocity. DOMAIN61 belongs to rejected Newton proposals on the inclined catenary; it does not demonstrate physical reverse flow in an accepted solution. The earlier positive-velocity derivative-stencil failure was removed. Bhagwat SS3 reports failure to detect an admissible recovered primitive root in the proposed states; this is not a global proof of physical root nonexistence. Exact cells, angles, trial factors and codes are retained in failure_diagnostics.json.

| OSC configuration | Period final 10 [s] | Semi-amplitude [Pa] | Pressure amplitude change [%] | 13-channel means/excursions <=1% |
|---|---:|---:|---:|---|
| N51, dt=0.125, end=2400 s | 23.10302275 | 9700.027290 | 2.397741 | True/False |
| N51, dt=0.0625, end=2400 s | 23.26563921 | 13763.965412 | 0.000523 | True/True |
| N101, dt=0.125, end=2400 s | 23.26639684 | 12769.699846 | 0.021282 | True/True |
| N51, dt=0.125, end=3600 s | 23.10864055 | 9909.563379 | 0.006414 | True/True |

Thirteen-channel means and excursions use original piecewise-linear time statistics on the final ten pressure-peak intervals. First-five versus last-five changes are divided by max(abs(first),abs(last),1e-30 in the channel unit); absolute changes and scales are recorded. Exact zero/zero gives zero. Pressure, two liquid holdups, four phase mass flows and six inventories are checked independently. Temporal accommodation does not imply convergence with mesh or step refinement.

The cycle observations show accumulation/compression proxies, liquid production and inventory depletion, with positive gas base flow throughout the observed cycles. Local front penetration is not directly measured. Complete gas blockage is not required or imposed for OSC/SS3. Neither oscillations alone nor unit-test passage constitute experimental validation. No parameters were fitted to force slugging.

The runner changes requested output_interval_s from 0.1 to 0.0625 with the temporal refinement; actual late monitor spacing is 0.125 versus 0.0625 s. The same piecewise-linear time statistics and quadratic extremum extraction are used, but the final-cycle signals are not observed on identical sampling grids. This observational difference is retained as a limitation of the temporal comparison. It changes recording cadence, not EOS, physical inputs or the integration step controller.

Conservation uses the original cumulative inlet/outlet phase-flow integrals. summary.json records absolute kg defects, normalization masses and normalized defects, including the stable control. Physical inputs (EOS, geometry, phase flows and fluid data) agree across closures and independent refinements; source/binary hashes stayed unchanged. The extension repeats all compared steps before the original shortened final
step at 2400 s; that last step is excluded from the bit-for-bit comparison.
The original 2400 s result remains preserved separately.

The independently rerun Bhagwat monitor files are compared byte-for-byte against the earlier 7321e67 cohort in bhagwat_regression_against_7321e67.json. This regression check does not merge the two physical-source cohorts.

The SS3 source-period conflict remains documented: dissertation 26 s; Baliño 2010, Table 3/p.653/case 11, 76 s. The present incomplete SS3 trajectories cannot resolve that discrepancy.

Full numerical evidence archive: `C:/Users/leonardoab/Documents/rotordin_frontend/artifacts/primsf-continuation-20261118/PRIMSF_Comparison_490936a_20261118.zip`
SHA256: `f90701bced6d634c971358dc4d9ceb09962907350cbbae69641a81eba31756c6` (46163553 bytes).
archive_audit.json records every CRC, all internal SHA256 checks, the ELF/PE signature scan and the allowed Matplotlib PDF figure list. No executable or bibliographic PDF is included.
