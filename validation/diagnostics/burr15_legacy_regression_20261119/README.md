# Real legacy/default regression — 2026-11-19

Both frozen versions were compiled and executed in separate detached worktrees:
`ee5810e` (legacy) and `8a85e37` (explicit fixed pipe policy). GNU Fortran 13.3.0,
release Makefile flags, `make -j1 FC=gfortran LIBS='-llapack -lblas'`; full commit,
source and binary SHA-256 values are in `summary.json`. Four solver executions
were sequential, with nice=15 and OMP/OpenBLAS/MKL thread counts fixed at one.
No campaign binary or source file was changed.

The common D_Stable and D_OSC benchmark decks came from ee5810e, using 51+51 cells,
the original descending pipe/catenary geometry, EOS, flows, Bhagwat riser and
MATCHED junction. Only duration was shortened identically to 1 s. Within each
pair the input text differs solely in out_dir, verified by a canonical input
SHA. Neither pipe_regime_mode nor closure_pipe is supplied: this exercises the
actual default policy. All four runs reached 1 s with status 0.

| Case | Complete monitor.csv | Complete fields.csv | max riser absolute Δul (m/s) | max riser absolute Δug (m/s) |
|---|---|---|---:|---:|
| D_Stable | byte-identical | byte-identical | 0.000111838 | 0.14066523 |
| D_OSC | byte-identical | byte-identical | 0.000397186 | 0.29012416 |

The ul/ug differences above belong only to the final legacy results file. The
new writer recovers its riser velocities with the selected closure ID3, whereas
the old writer omitted the ID and silently used Bendiksen. All pipe columns,
positions, pressure and alpha in that final file agree exactly at its printed
precision. The separate CLI suite compares final output to operator fields.
Byte equality here covers only these short trajectories; it is not long-time
or experimental validation.

The first harness attempt preserved a completed legacy D_Stable run, but its
long out_dir caused the legacy fixed-length filename buffer to truncate `.dat`
to `.d`. The harness stopped while locating that file. The final comparison used
new shorter output directories for both versions; the original attempt/log is
retained and no solver change was made to resolve this harness path issue.

`run_comparison.py` is the executed harness; it builds isolated worktrees and
asserts absent output directories to preserve prior evidence. Reproduction
requires new worktree/output labels rather than deleting or overwriting the
recorded data. Input decks, metadata, run status and solver logs are in the external evidence directory;
raw monitor/fields paths and SHA values are recorded in the JSON.

The compiler logs are preserved outside Git in `PRIMS-F-burr15-results-20261119/operator_evidence/legacy_build_logs`.

Additional versions, logs and raw records are preserved outside Git in
`PRIMS-F-burr15-results-20261119/operator_evidence/burr15_legacy_regression_20261119`.
