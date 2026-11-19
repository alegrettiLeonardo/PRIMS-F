from pathlib import Path
import csv
import json
import shutil

repo = Path('/home/leonardoab/Documentos/PRIMS-F-burr15-20261119')
report = repo / 'validation/diagnostics/burr15_rejection_replay_20261119'
provenance = json.loads((report / 'provenance_v2.json').read_text())


def values(lines, prefix):
    return [line[len(prefix):].split() for line in lines if line.startswith(prefix)]


def profile_data(lines, profile):
    result = {}
    prefixes = {
        'junction_SI': ('junction pstar_Pa,jl_m_s,jg_m_s profile ', ['p_Pa', 'jl_m_s', 'jg_m_s']),
        'junction_gas_terms_kg_s': ('junction gas convective,diffusive,net kg_s profile ', ['convective', 'diffusive', 'net']),
        'neighbors_SI': ('neighbor velocities ulL,ugL,ulR,ugR,alphaL,alphaR profile ', ['ulL_m_s', 'ugL_m_s', 'ulR_m_s', 'ugR_m_s', 'alphaL', 'alphaR']),
    }
    for key, (prefix, names) in prefixes.items():
        row = next(row for row in values(lines, prefix) if int(row[0]) == profile)
        result[key] = dict(zip(names, map(float, row[1:])))
    row = next(row for row in values(lines, 'riser trace root profile,status,nroots,alpha ') if int(row[0]) == profile)
    result['riser_trace_root'] = dict(status=int(row[1]), nroots=int(row[2]), alpha=float(row[3]) if int(row[1]) == 0 else None)
    return result


summary = {
    'schema': 'PRIMS-F.burr15_rejection_replay.v1',
    'frozen_head': provenance['frozen_head'],
    'campaign_binary_sha256': provenance['frozen_binary_sha256'],
    'source_and_binary_match_campaign': True,
    'snapshot_semantics': 'Uold is last accepted; Utrial is the final Newton linearization, not a saved rejected proposal',
    'proposal_norm_semantics': 'null when residual assembly failed; Fortran huge sentinel is not a measured residual',
    'runs': {},
}
for case in ('D_SS1_fixed', 'D_SS3_fixed'):
    item = {'snapshot_sha256': provenance['runs'][case]['snapshot_sha256'], 'modes': {}}
    for mode in ('F', 'T'):
        directory = report / (case + '_' + mode + '_v2')
        lines = (directory / 'probe.log').read_text().splitlines()
        snapshot = values(lines, 'snapshot time_s,dt_s,norm,status ')[0]
        preparation = int(values(lines, 'prepare_dissipation status ')[0][0])
        assert preparation == 0
        residual_difference = float(values(lines, 'residual difference from snapshot ')[0][0])
        assert residual_difference == 0
        assert not values(lines, 'point failure profile,cell,status ')
        assert not values(lines, 'face failure profile,leftcell,status ')
        assert all(int(row[1]) == 0 for row in values(lines, 'inlet boundary profile,status '))
        assert all(int(row[1]) == 0 for row in values(lines, 'outlet boundary profile,status '))
        proposals = []
        point_failures = values(lines, 'first point failure cell,status ')
        face_failures = values(lines, 'first face failure leftcell,status ')
        for index, row in enumerate(values(lines, 'proposal k,lambda,status,norm,Armijo ')):
            status = int(row[2])
            proposal = {
                'k': int(row[0]), 'lambda': float(row[1]), 'assembly_status': status,
                'residual_norm': float(row[3]) if status == 0 else None,
                'accepted_by_Armijo': row[4] == 'T',
                'first_point_failure_cell': int(point_failures[index][0]),
                'first_point_failure_status': int(point_failures[index][1]),
                'junction': profile_data(lines, 10 + index),
            }
            if face_failures:
                proposal['first_face_failure_leftcell'] = int(face_failures[index][0])
                proposal['first_face_failure_status'] = int(face_failures[index][1])
            proposals.append(proposal)
        assert len(proposals) == 18 and not any(p['accepted_by_Armijo'] for p in proposals)
        mode_summary = {
            'time_s': float(snapshot[0]), 'dt_s': float(snapshot[1]),
            'saved_residual_norm': float(snapshot[2]), 'terminal_status': int(snapshot[3]),
            'lagged_dissipation_status': preparation,
            'snapshot_residual_max_absolute_difference': residual_difference,
            'all_old_and_trial_points_faces_boundaries_ok': True,
            'old_junction': profile_data(lines, 1), 'linearization_junction': profile_data(lines, 2),
            'proposals': proposals,
        }
        if case == 'D_SS3_fixed':
            rows = list(csv.DictReader((directory / 'point_states.csv').open()))
            mode_summary['cell69_old_and_linearization'] = [
                {key: (int(value) if key in ('cell', 'profile', 'status') else float(value)) for key, value in row.items()}
                for row in rows if int(row['cell']) == 69
            ]
        item['modes'][mode] = mode_summary
    summary['runs'][case] = item
(report / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
shutil.copyfile(Path(__file__), report / 'summarize_replay.py')
shutil.copyfile(Path(__file__).with_name('build_burr15_rejection_replay.py'), report / 'build_replay.py')
(report / 'README.md').write_text('''# New Burr15-v1 rejection replay — frozen 8a85e37

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
''')
print('Recorded new SS1/SS3 rejection diagnosis with exact residual replay.')
