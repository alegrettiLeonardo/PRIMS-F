# P0/P1 qualification, 2026-11-19

This report is being completed against the new isolated branch. Historical data
retain their original operator identifiers; 05ad799 and 490936a are not combined.
Opening a draft MR is authorized; merging is not.

## Protocol fixed before the new transients

P0: inlet mass targets and accepted fluxes agree to 1e-12 relative for nonzero
targets; zero targets use absolute error and undefined relative error (NaN).
Outlet pressure is the pressure in the actual boundary flux after subtracting
its own kinetic term. Liquid pressure roundoff is bounded by the conditioning
of the affine EOS (64 epsilon rho c_l^2); no pressure is adjusted to pass.
Six phase/section inventory balances must satisfy |error|/M(0) <= 1e-6;
an exactly absent phase uses its absolute balance, no artificial denominator.
Phase fluxes on both sides of the internal junction must be identical.

Hydrostatics: one fixed catenary geometry, independently refined meshes,
G=0, 4.905, 9.81, 19.62 m/s2, exact absent gas inventory. Compare conditioned
roundoff residuals, velocities and accepted balances. Gas reduction is tested
independently. Full boundary Jacobians use independently perturbed residuals
with the ghost/face state recomputed and a step-size convergence curve.

SS1/SS3: first run to 110 s, a wall-time cap distinguishes truncation from
physical failure; preserve the last accepted checkpoint and rejected proposals.
Only sustained admissible advancement can remove the historical block. A
successful replay does not qualify ten cycles or a physical regime.

OSC: start with 2400 s, tol_res=tol_inc=1e-10, N51+51, dt_max=0.125 s.
The last ten complete peak-to-peak cycles must show <=0.5% variation in period,
amplitude and mean before use in refinement. Numerical convergence requires
independent h/dt changes and <=1% differences in period and semi-amplitude;
settling alone is not convergence. BAROTROPIC and ISOTHERMAL campaigns are
separate constitutive models; the latter uses c^2=RT, with no gamma.

P1 prototypes remain outside production. A source-compatible exterior trace
does not balance the entire mixed interior FV operator. The reduced operator
supports only globally absent inventory, not a moving mixed/pure interface.
Missing external reservoir data and unsupported branches return nonzero status.
No severe-slugging or experimental validation follows from these tests alone.
