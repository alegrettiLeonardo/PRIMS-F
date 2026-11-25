#!/bin/bash
# Global assembled-Jacobian check at accepted saved states (fields.csv) of the
# P1 campaign. Usage: run_global_jacobian.sh BINARY RUNS_DIR EVIDENCE_DIR
B=$1; R=$2; E=$3; D=$(dirname "$0"); mkdir -p "$E"
while read run t cid; do
  python3 "$D/prepare_state.py" "$R/$run" "$t" "$E/gj_${run}_t$t" >/dev/null || { echo "$run prepare FAIL"; continue; }
  "$B" "$E/gj_${run}_t$t.params" "$E/gj_${run}_t$t.state" "$cid" 2>&1 | grep -v 'runtime warning\|^At line' > "$E/gj_${run}_t$t.log"
  echo "$run t=$t closure=$cid $(grep JUNCTION.status "$E/gj_${run}_t$t.log" | cut -c1-60) | $(grep SUMMARY "$E/gj_${run}_t$t.log")"
done <<LIST
D_OSC_C5_BAROTROPIC_N51_dt0.5 60 5
D_OSC_C5_BAROTROPIC_N51_dt0.5 135 5
D_Stable_C5_BAROTROPIC_N51_dt0.5 600 5
D_SS1_C5_BAROTROPIC_N51_dt0.5 30 5
D_SS3_C5_BAROTROPIC_N51_dt0.5 30 5
D_OSC_C3_BAROTROPIC_N51_dt0.5 60 3
LIST
