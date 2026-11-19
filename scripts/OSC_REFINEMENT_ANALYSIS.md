# Cycle observations and conservation

These utilities post-process saved runs; they do not execute or modify the solver.

From the repository root, run the 17 regression cases with:

```sh
PYTHONPATH=scripts python3 -m unittest discover -s tests -p test_apparent_order.py
PYTHONPATH=scripts python3 -m unittest discover -s tests -p test_cycle_conservation.py
PYTHONPATH=scripts python3 -m unittest discover -s tests -p test_fixed_trace_holdup.py
```

`apparent_order.py` estimates observed order for three decreasing scales, including unequal mesh ratios. It keeps nonmonotonic differences, nonpositive orders, and differences below the declared resolution floor explicit. The caller must establish common physics, complete accommodated cycle windows, and effective time steps. An observed order does not establish an asymptotic range, resolution of the other axis, or experimental validity.

`inspect_actual_mesh.py --plan PLAN --output NEW.json` checks the initial saved FV centers against the declared uniform mesh and reconstructs faces; its scale is `L/N`. It refuses nonuniform declared meshes.

```sh
python3 scripts/analyze_cycle_conservation.py --run-dir RUN \
  --cycle-summary CYCLES.json --mass-summary P2_V2.json --output NEW_BALANCES.json
python3 scripts/analyze_fixed_trace_holdup.py --run-dir RUN \
  --cycle-summary CYCLES.json --output NEW_TRACE.json
```

The cycle summary must identify the run, preserve hashes of its monitor, manifest, metadata and status, and provide eleven edges for ten complete intervals. The conservation command additionally requires the accepted-ledger audit schema `PRIMS-F.accepted_conservation.v2`, with matching hashes. Its output keeps arithmetic-contract approval distinct from six physical phase/segment balances. It integrates signed accepted transport, uses absolute transport only as a normalization, and never subtracts Newton residuals from the physical error. Edge inventories and cumulative integrals are linearly interpolated; every accepted interior row participates in the maximum error. The inventory and processed-mass criteria remain `1e-6` without an added absolute allowance.

The supplementary holdup is `1-alpha_riser_trace` at the fixed junction. It uses the same ten declared cycle windows and 256 phase samples without changing the thirteen primary accommodation gates. It differs from `holdup_liquid_base`, which is measured at the first riser cell center and moves when the mesh changes. `cycle_waveform.py` retains fixed-phase error separately from a diagnostic lag; positive lag compares `fine(phi+lag)` with `coarse(phi)`. A constant waveform has no defined lag.

The campaign executed preserved external copies of these tools. Repository integration keeps the order, cycle-balance and mesh code byte-identical. The holdup wrapper now imports repository-local functions, reads the monitor directly after verifying the cycle-summary hashes, and records its current dependencies. The three interpolation/statistics functions in `cycle_waveform.py` were extracted unchanged from the executed analyzer; `validation_evidence.py` supplies the identical interpolation and monitor-loading implementations. Integrated tests and source hashes are recorded separately from the campaign's original tool manifests, so later repository edits cannot retroactively relabel its evidence.
