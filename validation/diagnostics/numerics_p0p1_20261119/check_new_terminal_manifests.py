"""Read-only terminal audit; no solver execution or data repair."""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import math
import re
import shutil

repo = Path('/home/leonardoab/Documentos/PRIMS-F-numerics-20261119')
campaign = repo.parent / 'PRIMS-F-numerics-results-20261119/v1_2400'
previous = repo.parent / 'PRIMS-F-burr15-results-20261119/v1/runs'
report = repo / 'validation/diagnostics/numerics_p0p1_20261119'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
plan = json.loads((campaign / 'plan.json').read_text())
seal = json.loads((campaign / 'seal.json').read_text())
assert sha(campaign / 'plan.json') == seal['plan_sha256']
assert sha(Path(seal['build_path'])) == seal['build_sha256']
assert sha(Path(seal['binary_path'])) == seal['binary_sha256']
frozen = Path(plan['root'])
assert all(sha(frozen / name) == digest for name, digest in seal['source_sha256'].items())
source = (frozen / 'scr/module_simulate_pipe_riser.f90').read_text()
allowed_stages = set(re.findall(r"call failure\([^\n]*?'([^'\n]+)'", source)) | {'NONE'}


def read_key_values(path):
    return dict(line.split('=', 1) for line in path.read_text().splitlines() if '=' in line)


summary = {
    'schema': 'PRIMS-F.numerics_terminal_manifest_check.v1',
    'frozen_head': seal['git_head'], 'binary_sha256': seal['binary_sha256'],
    'plan_sha256': seal['plan_sha256'], 'seal_sha256': sha(campaign / 'seal.json'),
    'all_source_hashes_verified': len(seal['source_sha256']),
    'plan_state_literal': plan['state'],
    'seal_verified': True, 'allowed_stage_enum_from_source': sorted(allowed_stages),
    'runs': {},
}
for name in ('D_SS1', 'D_SS3'):
    run = campaign / 'runs' / name
    old = previous / (name + '_fixed')
    manifest = json.loads((run / 'manifest.json').read_text())
    job = next(job for job in plan['jobs'] if job['name'] == name)
    assert manifest['git_head'] == seal['git_head']
    assert manifest['binary_sha256'] == seal['binary_sha256']
    assert manifest['source_sha256'] == seal['source_sha256']
    assert manifest['plan_sha256'] == seal['plan_sha256']
    assert manifest['seal_sha256'] == sha(campaign / 'seal.json')
    assert manifest['namelist_sha256'] == job['deck_sha256'] == sha(Path(job['deck']))
    assert manifest['binary_unchanged'] and manifest['sources_unchanged']
    status = {key: value.strip() for key, value in read_key_values(run / 'run_status.txt').items()}
    sidecar = {key: value.strip() for key, value in read_key_values(run / 'last_rejected_iteration.txt').items()}
    assert int(status['status']) == 6000 and manifest['returncode'] == 3
    assert status['termination_reason'] == sidecar['reason'] == 'NEWTON_GLOBALIZATION_FAILURE'
    with (run / 'rejected_proposals.csv').open() as stream:
        proposals = list(csv.DictReader(stream))
    for row in proposals:
        assert row['stage'] in allowed_stages
        index = int(row['extended_index'])
        if row['stage'] == 'CELL':
            assert 3 <= index <= 104
        if row['stage'] == 'FACE':
            assert 3 <= index <= 103
        if row['stage'] == 'NONE':
            assert index == 0 and int(row['assembly_status']) == 0
        assert row['full_jacobian'] in {'T', 'F'}
        assert 1 <= int(row['iteration']) <= 25
        assert 1 <= int(row['backtrack']) <= 18
        assert 0 < float(row['lambda']) <= 1
        norm = float(row['proposal_norm'])
        if int(row['assembly_status']) != 0:
            assert math.isnan(norm)
        else:
            assert math.isfinite(norm)
    last = proposals[-1]
    assert int(last['iteration']) == int(sidecar['newton_iteration'])
    assert int(last['backtrack']) == int(sidecar['backtrack']) == 18
    assert float(last['lambda']) == float(sidecar['lambda']) == 2**-17
    assert int(last['assembly_status']) == int(sidecar['proposal_status'])
    assert int(last['extended_index']) == int(sidecar['proposal_index'])
    assert last['stage'] == sidecar['proposal_stage']
    assert last['full_jacobian'] == sidecar['jac_dissipation_fd'] == 'T'
    assert math.isnan(float(sidecar['proposal_norm']))
    terminal = proposals[-18:]
    assert [int(row['backtrack']) for row in terminal] == list(range(1, 19))
    assert all(row['time_s'] == last['time_s'] and row['attempt_dt_s'] == last['attempt_dt_s']
               and row['iteration'] == last['iteration'] for row in terminal)
    assert all(row['assembly_status'] == last['assembly_status'] and row['stage'] == last['stage']
               and row['extended_index'] == last['extended_index'] for row in terminal)
    old_status = read_key_values(old / 'run_status.txt')
    assert float(status['reached_time_s']) == float(old_status['reached_time_s'])
    assert int(status['accepted_steps']) == int(old_status['accepted_steps'])
    assert int(status['rejected_steps']) == int(old_status['rejected_steps'])
    snapshot_equal = (run / 'last_rejected_state.dat').read_bytes() == (old / 'last_rejected_state.dat').read_bytes()
    assert snapshot_equal
    assert float(last['time_s']) == float(status['reached_time_s'])
    summary['runs'][name] = {
        'manifest_sha256': sha(run / 'manifest.json'), 'snapshot_sha256': sha(run / 'last_rejected_state.dat'),
        'snapshot_byte_equal_to_burr15_v1': snapshot_equal,
        'run_status_sha256': sha(run / 'run_status.txt'),
        'rejected_proposals_sha256': sha(run / 'rejected_proposals.csv'),
        'sidecar_sha256': sha(run / 'last_rejected_iteration.txt'),
        'reached_time_s': float(status['reached_time_s']),
        'accepted_steps': int(status['accepted_steps']), 'rejected_steps': int(status['rejected_steps']),
        'termination_reason': status['termination_reason'], 'terminal_status': int(status['status']),
        'iteration': int(sidecar['newton_iteration']), 'backtrack': int(sidecar['backtrack']),
        'lambda': float(sidecar['lambda']), 'full_jacobian': sidecar['jac_dissipation_fd'],
        'terminal_proposal_status': int(sidecar['proposal_status']),
        'terminal_proposal_index': int(sidecar['proposal_index']),
        'terminal_proposal_stage': sidecar['proposal_stage'], 'terminal_proposal_norm': None,
        'terminal_18_consistent': True, 'rejected_proposals_rows': len(proposals),
        'all_enum_and_index_checks_passed': True,
        'stage_status_counts': {str(key): value for key, value in Counter(
            (row['stage'], int(row['assembly_status'])) for row in proposals).items()},
        'raw_run': str(run),
    }
(report / 'terminal_69813ee_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
shutil.copyfile(Path(__file__), report / Path(__file__).name)
print(json.dumps({'head': summary['frozen_head'], 'sources_checked': summary['all_source_hashes_verified'],
                  'runs': {name: {key: run[key] for key in ['snapshot_sha256', 'snapshot_byte_equal_to_burr15_v1',
                  'terminal_proposal_status', 'terminal_proposal_index', 'terminal_proposal_stage', 'iteration',
                  'lambda', 'rejected_proposals_rows']} for name, run in summary['runs'].items()}}, indent=2))
