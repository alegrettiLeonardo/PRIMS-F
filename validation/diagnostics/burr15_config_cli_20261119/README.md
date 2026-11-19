# Fixed pipe configuration and CLI checks — 2026-11-19

Implemented `pipe_regime_mode` independently of the scenario `regime` label.
The default remains `FIXED_STRATIFIED`, with closure1. `FIXED_DRIFT` requires an
explicit `closure_pipe` ID2..7. `AUTO_BURR15` uses the stratified operator while
diagnostics are collected, and refuses a drift pipe closure. This configuration
does not select or manufacture a transient flow-regime transition.

Main passes the resolved pipe closure to both configuration and steady
initialization. The optional `closure_id_pipe` argument was appended to
`fun_dpds`, `steady_state_ndim_simp`, and `EstadoEstacionario_ndim_simp`, preserving
the positional contracts of older calls. Pipe section selection uses section
identity and is valid at beta=0. Pressure integration and reconstructed centre
states use the same pipe closure. No friction law or material parameter changed.

The final legacy `results_*.dat` reconstruction now uses explicit closure IDs
by section and checks recovery status. Previously its riser always defaulted
to Bendiksen even when the simulated closure was Bhagwat/Chexal. Thus ul/ug
columns in that post-processing file can change for old decks. Conservative
trajectories and the operator-based monitor/fields are not changed by this
output correction. Five CLI cases verify agreement between final results and
fields at the legacy output precision.

Executed checks, using a separate O0/checks build copied into
`tmp/pipe_config_check` so no integrated build objects were touched:

- Existing steady suite: **398/398 PASS**, including calls omitting the new
  optional argument.
- CLI: **18/18 PASS**. Five tiny 0.002s runs include default, explicit
  stratified, diagnostic AUTO, driftID2 at beta=0 and beta=2; six invalid
  configurations are refused; horizontal countercurrent Chexal fails separately
  during initialization with domain61. Default/stratified/AUTO physical monitor
  records are identical. All five final results match operator fields.

The first CLI version passed13 checks. Adding the five output comparisons
initially exposed the comparator's assumption of a unique final field frame;
the ee5810e writer repeats that frame. The comparator now retains the last
record at each cell position. No solver change or output-data repair was used.
The initial comparator-failure summary is retained by version.

Reproduce after building the solver:

```sh
python3 tests/test_pipe_regime_cli.py --binary bin/prims_f --output tmp/pipe_regime_cli_new
```

These are implementation/contract checks on tiny cases, not experimental
validation. Source hashes, executed logs and compiler provenance are retained.


The final v4 test separates deck, run-status, monitor and output-comparison
helpers into module functions, expands the approval conditions, and uses
explicit imports. The unchanged executable was tested again in a separate
directory: **18/18 PASS**. v3 summary and hashes are retained. The actual
ee5810e-to-8a85e37 short default regression is recorded separately in
`../burr15_legacy_regression_20261119/README.md`.

Additional versions, logs and raw records are preserved outside Git in
`PRIMS-F-burr15-results-20261119/operator_evidence/burr15_config_cli_20261119`.
