"""Run one predeclared smoke from an existing, verified frozen build."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import re
import resource
import shutil
import subprocess
import time


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def last_complete_row(path):
    if not path.exists():
        return None
    last = None
    with path.open() as stream:
        complete_lines = (line for line in stream if line.endswith('\n'))
        for row in csv.DictReader(complete_lines):
            if None not in row and None not in row.values():
                last = row
    return last


def parse_status(path):
    if not path.exists():
        return None
    result = {}
    for line in path.read_text().splitlines():
        if '=' in line:
            key, value = line.split('=', 1)
            result[key.strip()] = value.strip()
    return result


def hash_values(value):
    """Collect SHA strings; exact binary digest must occur in supplied manifest."""
    if isinstance(value, dict):
        return [item for child in value.values() for item in hash_values(child)]
    if isinstance(value, list):
        return [item for child in value for item in hash_values(child)]
    if isinstance(value, str) and re.fullmatch('[0-9a-fA-F]{64}', value):
        return [value.lower()]
    return []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--buildmanifest', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--case', choices=('D_SS1', 'D_SS3'), required=True)
    parser.add_argument('--timeout', type=float, default=900.0)
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--deck-root', type=Path, default=Path(
        '/home/leonardoab/Documentos/PRIMS-F-numerics-results-20261119/v1_2400/decks'))
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error('--timeout must be positive')

    binary = args.binary.resolve(strict=True)
    build_file = args.buildmanifest.resolve(strict=True)
    source_root = args.source_root.resolve(strict=True)
    build = json.loads(build_file.read_text())
    if build.get('qualified_build') is not True or build.get('returncode') != 0:
        raise SystemExit('The supplied build manifest is not qualified.')
    binary_hash = sha256(binary)
    declared_hashes = set(hash_values(build))
    if binary_hash not in declared_hashes:
        raise SystemExit('Binary SHA256 absent from the supplied frozen build manifest.')
    sources = {str(path.relative_to(source_root)): sha256(path)
               for path in sorted((source_root / 'scr').glob('*.f90'))}
    if not sources or sources != build.get('source_sha256'):
        raise SystemExit('Current Fortran source hashes differ from the frozen build manifest.')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source_root, text=True).strip()
    if commit != build.get('git_head') or commit != build.get('git_head_after'):
        raise SystemExit('Current source commit differs from the frozen build commit.')
    dirty = subprocess.check_output(['git', 'diff', '--name-only', 'HEAD', '--', 'scr'],
                                    cwd=source_root, text=True).strip()
    if dirty:
        raise SystemExit('Frozen source tree has local source changes: ' + dirty)
    if args.validate_only:
        print(json.dumps({'qualified': True, 'source_commit': commit,
                          'binary_sha256': binary_hash, 'case': args.case,
                          'timeout_wall_s': args.timeout, 'sources': len(sources)}))
        return

    duration = 100 if args.case == 'D_SS1' else 110
    folder = args.output_root.resolve() / args.case
    folder.mkdir(parents=True, exist_ok=False)
    output = folder / 'output'
    output.mkdir()
    local_binary = folder / 'prims_f'
    shutil.copy2(binary, local_binary)
    shutil.copy2(build_file, folder / 'build_manifest.json')
    source_deck = args.deck_root / (args.case + '.nml')
    deck_text = source_deck.read_text()
    deck_text, count_time = re.subn(r'time_dim\s*=\s*[0-9.eEdD+-]+',
                                   f'time_dim={duration}.0', deck_text)
    deck_text, count_out = re.subn(r"out_dir\s*=\s*'[^']*'",
                                  f"out_dir='{output}'", deck_text)
    if count_time != 1 or count_out != 1 or 'newton_globalization' in deck_text:
        raise SystemExit('Unexpected source deck structure; refusing ambiguous modifications.')
    deck_text, count_tolerance = re.subn(r'\btol_res\s*=\s*[0-9.eEdD+-]+',
                                        'tol_res=1e-10', deck_text)
    if count_tolerance != 1:
        raise SystemExit('Expected exactly one residual tolerance.')
    deck_text, count_end = re.subn(r'(?m)^\s*/\s*$',
                                 "  newton_globalization='ACTIVE_SET'\n/", deck_text)
    if count_end != 1:
        raise SystemExit('Expected exactly one namelist terminator.')
    deck = folder / (args.case + '.nml')
    deck.write_text(deck_text.rstrip() + '\n')
    environment = os.environ.copy()
    thread_environment = {name: '1' for name in
                          ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}
    environment.update(thread_environment)
    command = [str(local_binary), str(deck)]
    manifest = {
        'schema': 'primsf_frozen_active_smoke_tolerance_v1', 'case': args.case,
        'requested_time_s': duration, 'timeout_wall_s': args.timeout,
        'source_commit': commit, 'source_root': str(source_root), 'source_hashes': sources,
        'binary_origin': str(binary), 'binary_sha256': binary_hash,
        'build_manifest_origin': str(build_file), 'build_manifest_sha256': sha256(build_file),
        'source_deck_sha256': sha256(source_deck), 'deck_sha256': sha256(deck),
        'deck_changes': ['time_dim', 'out_dir', "newton_globalization='ACTIVE_SET'", 'tol_res=1e-10'],
        'tolerance_sensitivity': {'tol_res': 1e-10, 'tol_inc_unchanged': 1e-8, 'reference_tol_res': 1e-8},
        'helper_sha256': sha256(Path(__file__)), 'command': command,
        'environment': thread_environment, 'start_utc': utc_now(), 'state': 'RUNNING',
        'harness_pid': os.getpid(), 'output_directory': str(output),
    }
    manifest_file = folder / 'manifest.json'
    start = time.monotonic()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with (folder / 'solver.log').open('w') as log:
        process = subprocess.Popen(command, cwd=source_root, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT)
        manifest['solver_pid'] = process.pid
        atomic_json(manifest_file, manifest)
        print(json.dumps({'case': args.case, 'pid': process.pid,
                          'manifest': str(manifest_file), 'started': manifest['start_utc']}), flush=True)
        timed_out = False
        try:
            returncode = process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.kill()
            returncode = process.wait()
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    status = parse_status(output / 'run_status.txt')
    reason = 'EXTERNAL_TIMEOUT' if timed_out else 'PROCESS_EXIT'
    if not timed_out and status is not None:
        reason = status.get('termination_reason', reason)
    manifest.update({
        'state': 'FINISHED', 'finish_utc': utc_now(),
        'elapsed_wall_s': time.monotonic() - start,
        'child_cpu_seconds': usage_after.ru_utime + usage_after.ru_stime
                             - usage_before.ru_utime - usage_before.ru_stime,
        'returncode': returncode, 'timed_out': timed_out, 'termination_reason': reason,
        'run_status': status, 'last_complete_step': last_complete_row(output / 'steps.csv'),
        'last_complete_conservation': last_complete_row(output / 'conservation.csv'),
        'last_rejected_iteration': parse_status(output / 'last_rejected_iteration.txt'),
        'snapshot_role': ('last prior rejection, not timeout state' if timed_out
                          else 'last recorded rejection; consult termination_reason'),
        'progress_role': ('last fully flushed record; lower bound at timeout' if timed_out
                          else 'last complete recorded step'),
        'files_sha256': {path.name: sha256(path) for path in sorted(output.iterdir()) if path.is_file()},
    })
    atomic_json(manifest_file, manifest)
    print(json.dumps({'case': args.case, 'termination_reason': reason,
                      'wall_s': manifest['elapsed_wall_s'], 'cpu_s': manifest['child_cpu_seconds'],
                      'last_complete_step': manifest['last_complete_step']}), flush=True)


if __name__ == '__main__':
    main()
