# Thirteen-channel accommodation supplement

The main 7321e67 archive remains unchanged. This supplement measures first-five
versus last-five cycle means and excursions for the same 13 channels as the frozen
baseline: pressure, two liquid holdups, four phase mass flows and six inventories.
Original piecewise-linear monitor statistics are integrated on the ten pressure-peak
intervals. Relative changes use max(abs(first),abs(last),1e-30 in channel units);
absolute changes and normalization scales are retained, so zero/near-zero signals
can be reviewed. Exact zero/zero yields zero.

The explicit 1% checks test every channel separately. Passing temporal checks does
not establish mesh/time-step convergence, physical regime or experimental validity.
See summary.json for the worst channels and cycle_physics.json for all 13 signals.
The supplement archive, hashes and CRCs are recorded in archive_manifest.json;
it contains only generated figures, numerical tables and this analysis source.
