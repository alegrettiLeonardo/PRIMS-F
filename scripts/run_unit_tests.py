#!/usr/bin/env python3
"""Clean Fortran-2018 verification build, excluding main; no external Python packages.

Usage (from any working directory):
    FC=gfortran LIBS='-llapack -lblas' python3 scripts/run_unit_tests.py

Production code does not require procedure trampolines. CALLBACK_FLAGS and
LDFLAGS may be overridden for an independently qualified compiler.
All required tests must exist. Missing tests are failures, never silent skips.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone

SOURCES = [
    'mod_kinds', 'mod_rootfinding', 'mod_void_helpers', 'module_flow_regime',
    'module_config', 'module_lapack_small', 'module_chexal_1992', 'module_closure_law', 'module_find',
    'module_catenary', 'module_steady_state', 'module_jacobian', 'module_matrix_roe',
    'module_roe_solver', 'module_feed_contracts', 'module_boundary_trace', 'module_boundary_conditions', 'module_solver_utils',
    'module_time_step', 'module_active_set', 'module_single_phase', 'module_simulate_pipe_riser',
]
TESTS = ['test_chexal_velocity_reference', 'test_chexal_derivative_tail', 'test_chexal_1992_core', 'test_chexal_integration', 'test_geometry', 'test_closure_physics', 'test_primitives_numerics', 'test_steady',
         'test_boundaries_time', 'test_fv_operators', 'test_roe_geometry', 'test_chexal_junction', 'test_burr15', 'test_active_set', "test_active_replay", "test_feed_contracts", "test_partition_invariance", "test_single_phase", "test_source_boundary", "test_initial_inventory", "test_reduced_domain_edges", "test_junction_trace", "test_hydrostatic_outlet"]

def run_command(command: list[str], cwd: Path, log: Path, timeout: float) -> dict:
    start = time.perf_counter()
    result = {'command': command, 'log': str(log), 'status': 'FAIL'}
    try:
        process = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, timeout=timeout, check=False)
        output = process.stdout
        result.update(returncode=process.returncode,
                      status='PASS' if process.returncode == 0 else 'FAIL')
    except subprocess.TimeoutExpired as exc:
        raw = exc.stdout or ''
        output = raw.decode(errors='replace') if isinstance(raw, bytes) else raw
        output += '\nTIMEOUT: process terminated by test runner.\n'
        result.update(returncode=None, status='TIMEOUT')
    except OSError as exc:
        output = f'{type(exc).__name__}: {exc}\n'
        result.update(returncode=None, status='BLOCKED')
    log.write_text('$ ' + shlex.join(command) + '\n' + output, encoding='utf-8')
    result['elapsed_seconds'] = time.perf_counter() - start
    result['output'] = output
    return result

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, default=None)
    parser.add_argument('--timeout', type=float, default=180.0, help='Per-test timeout, seconds')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = (args.build_dir or root / 'tmp' / 'unit_tests').resolve()
    build.mkdir(parents=True, exist_ok=True)
    logs = build / 'logs'
    logs.mkdir(exist_ok=True)
    summary_path = build / 'summary.json'
    fc = shlex.split(os.environ.get('FC', 'gfortran'))
    libs = shlex.split(os.environ.get('LIBS', '-llapack -lblas'))
    callbacks = shlex.split(os.environ.get('CALLBACK_FLAGS', ''))
    link_flags = shlex.split(os.environ.get('LDFLAGS', '-fopenmp -Wl,-z,noexecstack'))
    flags = ['-std=f2018', '-O0', '-g3', '-Wall', '-Wextra',
             '-Wno-unused-dummy-argument', '-Wno-unused-parameter',
             '-fimplicit-none', '-ffree-line-length-none', '-fopenmp',
             '-fcheck=all', '-fbacktrace', '-ffpe-trap=invalid,zero,overflow',
             '-finit-real=snan', '-finit-integer=-999', *callbacks,
             '-J' + str(build), '-I' + str(build), '-I' + str(root / 'scr')]
    summary = {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'FAIL',
               'compiler': fc, 'flags': flags, 'link_flags': link_flags, 'libraries': libs,
               'build_dir': str(build), 'build': [], 'tests': [], 'source_sha256': {}}
    try:
        summary['git_head'] = subprocess.check_output(
            ['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True, timeout=10).strip()
    except (OSError, subprocess.SubprocessError):
        summary['git_head'] = None
    required = [(root / 'scr' / (name + '.f90')) for name in SOURCES]
    required += [(root / 'tests' / (name + '.f90')) for name in TESTS]
    replay_support = root / 'validation/diagnostics/qualification2_20261119/globalization_flux'
    required += [replay_support / 'offline_api.f90',
                 replay_support / 'fixtures/D_SS1_active_terminal/last_rejected_state.dat',
                 replay_support / 'fixtures/D_SS1_active_terminal/probe_parameters.dat',
                 replay_support / 'fixtures/D_SS3_active_terminal/last_rejected_state.dat',
                 replay_support / 'fixtures/D_SS3_active_terminal/probe_parameters.dat']
    missing = [str(path.relative_to(root)) for path in required if not path.is_file()]
    if missing:
        summary['status'] = 'MISSING_REQUIRED_FILES'
        summary['missing'] = missing
        summary_path.write_text(json.dumps(summary, indent=2) + '\n')
        print('FAIL: required files absent: ' + ', '.join(missing), file=sys.stderr)
        return 2
    if not fc or args.timeout <= 0:
        summary['status'] = 'INVALID_RUNNER_CONFIGURATION'
        summary_path.write_text(json.dumps(summary, indent=2) + '\n')
        return 2
    for path in required:
        summary['source_sha256'][str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    # Remove only objects/modules owned by this runner; do not reuse stale builds.
    for name in SOURCES:
        for suffix in ('.o', '.mod'):
            (build / (name + suffix)).unlink(missing_ok=True)
    version = run_command([*fc, '--version'], root, logs / 'compiler.log', 30)
    summary['compiler_version'] = version.pop('output')
    summary['compiler_check'] = version
    failed = version['status'] != 'PASS'
    objects: list[str] = []
    if not failed:
        for name in SOURCES:
            obj = build / (name + '.o')
            command = [*fc, *flags, '-c', str(root / 'scr' / (name + '.f90')), '-o', str(obj)]
            entry = run_command(command, root, logs / ('compile_' + name + '.log'), 180)
            entry.pop('output')
            entry['module'] = name
            summary['build'].append(entry)
            if entry['status'] != 'PASS':
                failed = True
                print(f"FAIL compile {name}; see {entry['log']}")
                break
            objects.append(str(obj))
    if not failed:
        for name in TESTS:
            executable = build / name
            executable.unlink(missing_ok=True)
            command = [*fc, *flags, str(root / 'tests' / (name + '.f90')), *objects,
                       *link_flags, *libs, '-o', str(executable)]
            link = run_command(command, root, logs / ('link_' + name + '.log'), 180)
            link.pop('output')
            entry = {'name': name, 'link': link, 'status': link['status']}
            if link['status'] == 'PASS':
                run_args = [str(executable)]
                if name == 'test_active_replay':
                    run_args += [str(root), str(build)]
                run = run_command(run_args, build, logs / (name + '.log'), args.timeout)
                output = run.pop('output')
                matches = re.findall(r'PASS[^\n]*?:\s*(\d+)\s+assertions', output)
                matches += re.findall(r'PASS\s+(\d+)\s+[^\n]*\bchecks\b', output)
                run['assertions'] = sum(map(int, matches)) if matches else None
                if run['status'] == 'PASS' and not re.search(r'\bPASS\b', output):
                    run['status'] = 'MISSING_PASS_MARKER'
                entry['run'] = run
                entry['status'] = run['status']
            summary['tests'].append(entry)
            print(f"{entry['status']}: {name}")
            # Independent tests continue so every failure has its own evidence.
        failed = any(test['status'] != 'PASS' for test in summary['tests'])
    summary['sources_changed_during_run'] = [
        str(path.relative_to(root)) for path in required
        if hashlib.sha256(path.read_bytes()).hexdigest() != summary['source_sha256'][str(path.relative_to(root))]
    ]
    if summary['sources_changed_during_run']:
        failed = True
        print('FAIL: sources changed during verification; rerun against a stable source snapshot.')
    summary['status'] = 'PASS' if not failed and len(summary['tests']) == len(TESTS) else 'FAIL'
    summary['passed_tests'] = sum(test['status'] == 'PASS' for test in summary['tests'])
    summary['required_tests'] = len(TESTS)
    summary['assertions'] = sum(test.get('run', {}).get('assertions') or 0 for test in summary['tests'])
    summary_path.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(f"{summary['status']}: {summary['passed_tests']}/{len(TESTS)} tests; summary: {summary_path}")
    return 0 if summary['status'] == 'PASS' else 1

if __name__ == '__main__':
    raise SystemExit(main())
