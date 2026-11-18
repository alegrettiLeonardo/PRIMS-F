#!/usr/bin/env python3
"""Test the candidate patch in temporary copies; never modify scr/ or bin/."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import shlex
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
MODULES = ('mod_kinds.f90', 'mod_void_helpers.f90', 'module_closure_law.f90',
           'module_steady_state.f90')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(argv, cwd, log, stdin=None):
    result = subprocess.run(argv, cwd=cwd, input=stdin, text=True, capture_output=True)
    log.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}); see {log.name}')
    return result.stdout


def compile_run(tree, program, mode, optimization, compiler, output, label, stdin=None):
    build = tree / ('build_' + label)
    build.mkdir()
    target = build / 'test'
    argv = compiler + ['-std=f2018', '-' + optimization, '-g', '-ffree-line-length-none',
                       '-fcheck=all,no-array-temps', '-ffpe-trap=invalid,zero,overflow',
                       '-J' + str(build), '-I' + str(build)]
    argv += [str(tree / 'scr' / name) for name in MODULES] + [str(program), '-o', str(target)]
    command(argv, tree, output / (label + '_compile.log'))
    return command([str(target), *mode], tree, output / (label + '_run.log'), stdin)


def review_inputs():
    fixtures = json.loads((HERE / 'tail_fixtures.json').read_text())['cases']
    scales = [(1., 1.), (997.05, .53), (1e6, 1e-6), (1e-6, 1e6), (997.05, 1e6), (1e6, .53)]
    cases = [dict(row, wr=wr, wu=wu, tol=tol, kind='positive') for row in fixtures
             for wr, wu in scales for tol in (1e-6, 1e-12, 1e-14)]
    cases += [dict(fixtures[1], jg_m_s=str(jg), wr=997.05, wu=.53, tol=1e-12, kind='excluded')
              for jg in (0., -1e-12, -1e-300)]
    inputs = ''.join(' '.join(map(str, [i, c['jl_m_s'], c['jg_m_s'], c['theta_rad'],
                                      c['wr'], c['wu'], c['tol']])) + '\n'
                     for i, c in enumerate(cases, 1))
    return cases, inputs


def compare_review(cases, old, new):
    old_lines, lines = old.splitlines(), new.splitlines()
    if len(old_lines) != len(cases) or len(lines) != len(cases):
        raise RuntimeError('Independent review output count differs from inputs')
    rows = []
    for case, before, after in zip(cases, old_lines, lines):
        f = after.split()
        index, code, nroots, s0, s1, rs = map(int, f[:6])
        a, h0, h1, residual = map(float, f[6:])
        old_status = int(before.split()[1])
        row = dict(index=index, status=code, baseline_status=old_status)
        if case['kind'] == 'positive':
            err = abs(a / float(case['alpha']) - 1)
            rerr = abs(residual / (float(case['jg_m_s']) / case['wu']))
            accurate = (code == s0 == s1 == rs == 0 and nroots == 1 and math.isfinite(a)
                        and a == h0 == h1 and err < 5e-12 and rerr < 5e-12)
            if old_status == 0:
                category, passed = 'preserved_old_success', before == after
            else:
                category, passed = 'new_lower_tail_root', old_status == 21 and accurate
            row.update(relative_alpha_error=err, relative_closure_residual=rerr,
                       accuracy_criterion_met=accurate)
        else:
            expected = 21 if float(case['jg_m_s']) == 0 else 12
            category = 'preserved_exclusion'
            passed = code == s0 == s1 == expected and nroots == 0 and math.isnan(a) and before == after
        row.update(category=category, patch_contract_met=passed)
        rows.append(row)
    if not all(r['patch_contract_met'] for r in rows):
        raise RuntimeError('Independent review contract failed')
    return dict(cases=len(rows), passed=len(rows), results=rows,
                preexisting_accuracy_limit_cases=sum(r['category'] == 'preserved_old_success'
                                                     and not r.get('accuracy_criterion_met', True) for r in rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--output', type=Path, help='Log directory; default tmp/alpha_tail_candidate_tests in the repository')
    parser.add_argument('--check-only', action='store_true', help='Check package hashes and patch application in a temporary copy; do not compile')
    parser.add_argument('--optimization', nargs='+', choices=('O0', 'O2', 'O3'), default=['O0', 'O3'])
    parser.add_argument('--review', action='store_true', help='Also reproduce the independent 129-case O2 scale/tolerance review')
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = (args.output or repo / 'tmp/alpha_tail_candidate_tests').resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((HERE / 'candidate_manifest.json').read_text())
    for name, expected in manifest['files'].items():
        if sha(HERE / name) != expected:
            raise RuntimeError(f'Candidate payload checksum differs: {name}')
    before_hashes = {name: sha(repo / 'scr' / name) for name in MODULES}
    binary = repo / 'bin/prims_f'
    binary_before = sha(binary) if binary.is_file() else None
    current = before_hashes['module_steady_state.f90']
    old_hash, new_hash = manifest['baseline_steady_sha256'], manifest['candidate_steady_sha256']
    if current not in (old_hash, new_hash):
        raise RuntimeError('The steady source is neither the reviewed baseline nor the exact candidate; rebase/review explicitly')
    compiler = shlex.split(os.environ.get('FC', 'gfortran'))
    report = {'candidate_package_status': 'CANDIDATE_NOT_APPLIED',
              'detected_worktree_state': 'baseline' if current == old_hash else 'candidate_already_applied',
              'check_only': args.check_only, 'patch_sha256': sha(HERE / 'alpha_tail.patch'),
              'production_source_before': before_hashes, 'runs': []}
    with tempfile.TemporaryDirectory(prefix='prims_alpha_tail_') as temporary:
        trees = {}
        for variant in ('before', 'candidate'):
            tree = Path(temporary) / variant
            (tree / 'scr').mkdir(parents=True)
            for name in MODULES:
                shutil.copyfile(repo / 'scr' / name, tree / 'scr' / name)
            shutil.copyfile(HERE / 'tail_fixtures.dat', tree / 'tail_fixtures.dat')
            command(['git', 'init', '-q', str(tree)], tree, output / (variant + '_git_init.log'))
            reverse = variant == 'before' and current == new_hash
            apply = variant == 'candidate' and current == old_hash
            if reverse or apply:
                argv = ['git', 'apply'] + (['--reverse'] if reverse else [])
                command(argv + ['--check', str(HERE / 'alpha_tail.patch')], tree,
                        output / (variant + '_patch_check.log'))
                command(argv + [str(HERE / 'alpha_tail.patch')], tree,
                        output / (variant + '_patch_apply.log'))
            expected = old_hash if variant == 'before' else new_hash
            if sha(tree / 'scr/module_steady_state.f90') != expected:
                raise RuntimeError('Temporary patched source hash is not the reviewed candidate')
            trees[variant] = tree
        report['temporary_patch_hash_verified'] = True
        if not args.check_only:
            for opt in args.optimization:
                emitted = {}
                for variant, tree in trees.items():
                    label = f'{variant}_{opt}_old_paths'
                    emitted[variant] = compile_run(tree, HERE / 'test_alpha_tail_guard.f90',
                                                   ['emit'], opt, compiler, output, label)
                if emitted['before'] != emitted['candidate'] or len(emitted['before'].splitlines()) != 101:
                    raise RuntimeError(f'Old roots/statuses differ under {opt}')
                text = compile_run(trees['candidate'], HERE / 'test_alpha_tail_guard.f90',
                                   ['check'], opt, compiler, output, f'candidate_{opt}_guards')
                if 'PASS 55 alpha tail guard checks' not in text:
                    raise RuntimeError('Expected 55 tail checks were not reported')
                for variant, tree in trees.items():
                    text = compile_run(tree, repo / 'tests/test_steady.f90', [], opt, compiler,
                                       output, f'{variant}_{opt}_steady')
                    if 'PASS 398' not in text:
                        raise RuntimeError('Expected 398 steady checks were not reported')
                report['runs'].append({'optimization': opt, 'tail_guards_pass': 55,
                                       'steady_checks_before_pass': 398, 'steady_checks_candidate_pass': 398,
                                       'old_records_byte_identical': 101})
            if args.review:
                cases, inputs = review_inputs()
                results = {variant: compile_run(tree, HERE / 'probe_alpha_tail_review.f90', [], 'O2',
                                               compiler, output, variant + '_O2_independent_review', inputs)
                           for variant, tree in trees.items()}
                report['independent_review'] = compare_review(cases, results['before'], results['candidate'])
    report['production_sources_unchanged'] = all(sha(repo / 'scr' / name) == value for name, value in before_hashes.items())
    report['production_binary_unchanged'] = (sha(binary) if binary.is_file() else None) == binary_before
    if not report['production_sources_unchanged'] or not report['production_binary_unchanged']:
        raise RuntimeError('Production changed while checks were running')
    report['status'] = 'CHECKED_WITHOUT_COMPILATION' if args.check_only else 'PASS'
    (output / 'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ('production_source_before', 'independent_review')}, indent=2))


if __name__ == '__main__':
    main()
