# New Burr15-v1 rejection replay — frozen 8a85e37

This is new evidence from `v1/runs/D_SS1_fixed` and `D_SS3_fixed`, using their
own exact conservative snapshots and metadata. The probe links release objects
from detached `8a85e37319995b37795d2c798d6f865b863b1834`; source hashes and the
solver binary SHA match the campaign manifest. The binary SHA is
`6b3164d2a18e08e173c1e65e5058bfc024d648a94b64fa10c7bc8ab8a24d37e5`.
No solver/campaign source or binary was changed. Only offline probes ran,
sequentially with nice=15 and single-threaded BLAS/OpenMP.

Both snapshots reproduce their saved residual vectors with maximum absolute
difference **zero**. Preparation of lagged Roe dissipation, all 102 physical
cells, all internal faces, both boundaries and Jacobian assembly succeed at
both the last accepted state and saved Newton linearization. Both F/T Jacobian
settings were replayed explicitly; every one of their 18 proposed steps fails
residual assembly before an Armijo norm can be evaluated. The proposals use the
actual positivity-limited starting lambda (1 here), then 1/2 through 2^-17.

| Case | Last accepted time (s) | Saved norm | First proposal refusal |
|---|---:|---:|---|
| SS1 | 78.571788476739314 | 2.9764063309210322e-7 | MATCHED face 53/54, riser trace, face1012 / assembly3012 |
| SS3 | 97.131852338700099 | 5.0692576609245183e-8 | Cell69 primitive recovery14 / assembly1014 |

SS1: old/common-trace jg is +3.9866594701254724e-9 m/s; at the saved
linearization it is +2.0694867733244484e-13 m/s. The corresponding riser roots
are positive, approximately 4.47e-9 and 2.32e-13. The next full Newton proposal
has p*=281717.21888713469 Pa, jl*=0.52948669929489744 m/s and
jg*=-6.0549846226809329e-9 m/s. At lambda=2^-17, jg* is still
-3.2402369683065269e-14 m/s. All physical-cell recoveries succeed for these
proposals; scanning faces locates the first refusal at face53, the pipe/riser
junction. The Bhagwat superficial-flux trace returns12, wrapped as1012 then3012.
Convective, dissipative and net gas mass fluxes in kg/s are retained in the JSON
and logs. The negative numerical trace does not establish physical gas reversal.

SS3: every proposed lambda first fails recovery at extended cell69 (physical
cell67, riser cell16), theta=0.90123866722684609 rad. Its last accepted state has
alpha=0.63827166303094773, p=228440.19275256811 Pa, ul=2.9367850460519817e-9 m/s
and ug=0.62455971208025673 m/s. The saved linearization still has positive
ul=2.2028346473388749e-14 m/s. Recovery NO_ROOT4 is wrapped as point14 then1014;
this means the algorithm did not find an admissible root and is not a global
proof of physical nonexistence. The junction remains admissible, with positive
jg near0.226935 m/s. This failure is distinct from SS1's negative gas trace.

The actual terminal code in both runs is6000: line-search failure after the
fallback. This frozen version does not save Newton iteration or lambda in a
sidecar. Therefore the reported lambdas are explicit offline reconstructions
from the saved linearization, not claimed telemetry or a reconstructed complete
Newton history. `Uold` is accepted; `Utrial` is the saved linearization rather
than the rejected proposal. Failed-assembly huge sentinels become null norms in
JSON and are never reported as measured residual magnitudes.

The original v1 probe and logs are preserved. The v2 probe adds only an offline
face scan identifying SS1 face53. Full source/object/snapshot hashes and
parameters are in provenance files; `build/` contains ignored local executables.
The scripts and probe sources reproduce the diagnosis without clipping, gas
projection, root substitution or any change to the frozen operator.

The compact tracked result is `compact_summary.json`. Repeated raw snapshots,
logs and full provenance/18-proposal tables remain local and are preserved
in the ZIP identified by `archive_manifest.json`; CRC and every member SHA
were verified. No PDF, ELF/PE binary or Fortran object/module is included.
