# Frozen pre-derivative-fix comparison (2026-11-18)

`7321e6744932200176a0a6628b07e195d342071c` contains alpha-tail, augmented Roe geometry, Chexal1992 and coupled adapters. It precedes the near-zero velocity derivative correction. This complete historical cohort remains separate from the corrected-source comparison.

All cases use BAROTROPIC / MATCHED / LAGGED / ROE, original geometry, EOS and imposed phase inputs. Source and binary hashes stayed unchanged; physical input equality across closures and independent refinements was checked. No fitted coefficients or parameter tuning were used.

| Case | Closure | Cells per section | dt max [s] | End [s] | Process result | Solver status |
|---|---:|---:|---:|---:|---|---:|
| D_OSC | 3 | 51 | 0.125 | 2400.00000000 | COMPLETED | 0 |
| D_OSC | 5 | 51 | 0.125 | 93.10000408 | FAILED | 6000 |
| D_SS1 | 3 | 51 | 0.125 | 78.57178848 | FAILED | 6000 |
| D_SS1 | 5 | 51 | 0.125 | 52.79526549 | FAILED | 3061 |
| D_SS3 | 3 | 51 | 0.125 | 97.13185234 | FAILED | 6000 |
| D_SS3 | 5 | 51 | 0.125 | 42.18200940 | FAILED | 6000 |
| D_Stable | 3 | 51 | 0.125 | 2400.00000000 | COMPLETED | 0 |
| D_Stable | 5 | 51 | 0.125 | 2400.00000000 | COMPLETED | 0 |
| D_OSC | 3 | 51 | 0.0625 | 2400.00000000 | COMPLETED | 0 |
| D_OSC | 5 | 51 | 0.0625 | 77.36499358 | FAILED | 6000 |
| D_OSC | 3 | 51 | 0.125 | 3600.00000000 | COMPLETED | 0 |
| D_OSC | 3 | 101 | 0.125 | 2400.00000000 | COMPLETED | 0 |
| D_OSC | 5 | 101 | 0.125 | 94.24528560 | FAILED | 6000 |

Closure 3 is Bhagwat; closure 5 is Chexal1992 AUTO. Failed runs were not assigned periods from truncated substitute windows.

Failure diagnosis: SS1 reaches a negative numerical gas trace in the horizontal MATCHED mapping while cell states remain admissible (Bhagwat trace 12; Chexal 61). SS3 Bhagwat proposals return NO_ROOT in the local cocurrent primitive recovery;
this finite search is not a global proof that no admissible root exists. Chexal SS3 and OSC/refinements stop at a positive liquid velocity near 2.2586e-8 m/s because the old derivative stencil cannot remain on the admissible piece. Some larger OSC Newton proposals are actually negative on the inclined catenary and correctly fail 61. The b0d33cb same-snapshot replays, recorded separately, remove the numerical 53 failures while preserving domain refusals. No projection or artificial spectrum was used.

Temporal accommodation uses the existing declared 1% amplitude-change check over ten complete pressure cycles. Conservation uses original cumulative phase-flow integrals; both absolute kg errors and normalization-relative errors are recorded. These are numerical criteria, not experimental validation.

| OSC configuration | T final10 [s] | Mean semi-amplitude [Pa] | Final10 amplitude change [%] |
|---|---:|---:|---:|
| N51, dt=0.125, end=2400 s | 23.10302275 | 9700.027290 | 2.397741 |
| N51, dt=0.0625, end=2400 s | 23.26563921 | 13763.965412 | 0.000523 |
| N101, dt=0.125, end=2400 s | 23.26639684 | 12769.699846 | 0.021282 |
| N51, dt=0.125, end=3600 s | 23.10864055 | 9909.563379 | 0.006414 |

The N51/dt.125 extension reproduces the original monitor bit-for-bit through 19201 rows, excluding only the original final shortened step at 2400 s. Its amplitude settles after extension to 3600 s, while the N101 amplitude differs materially; temporal accommodation does not establish mesh convergence. See cycles.json for exact windows, threshold flags and period/amplitude variation.

Cycle-conditioned signals include pressure, liquid fractions, pipeline gas and riser liquid inventories, and both phase fluxes at base/outlet. Synchronous pressure/inventory increase supports accumulation/compression proxies; excess liquid production accompanies riser depletion; pressure decrease accompanies gas-inventory depletion. Gas flow at the base remains positive. These observations support a numerical oscillatory mechanism; local gas-front penetration was not measured directly and no universal complete-blockage criterion is imposed on OSC/SS3.

The stable control stays flat over 900–2400 s with both closures; differing equilibrium pressure/holdup is reported rather than tuned away. No failed transient or pressure oscillation is treated as experimental validation. SS3 reference-period disagreement remains: dissertation 26 s versus Baliño 2010 Table 3/p.653/case 11: 76 s.

Full series, fields, effective namelists, analysis figures and detailed logs are preserved locally without executables or bibliographic PDF content:
- Archive: `C:/Users/leonardoab/Documents/rotordin_frontend/artifacts/primsf-continuation-20261118/PRIMSF_Comparison_7321e67_20261118.zip`
- SHA256: `6b05b55e46bb8be86bd8f3f8245b3337809b860c421dc0fdc1c1254ef4bc6808`
- Bytes: 47591863
- CONTENTS_SHA256.json inside the archive records every included file hash. Compact metrics and reproducible analysis scripts are versioned here.
