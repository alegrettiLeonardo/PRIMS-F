from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

base = Path('/home/leonardoab/Documentos/PRIMS-F-numerics-20261119')
frozen = base.parent / 'PRIMS-F-burr15-regression-8a85e37-20261119'
report = base / 'validation/diagnostics/numerics_p0p1_20261119'
build = base / 'tmp/p0p1_probe'
report.mkdir(parents=True, exist_ok=True)
build.mkdir(parents=True, exist_ok=True)
src = Path(__file__).with_name('probe_interfaces.f90')
shutil.copyfile(src, report / src.name)
shutil.copyfile(Path(__file__), report / Path(__file__).name)
objects = sorted(p for p in (frozen / 'tmp/release').glob('*.o') if p.name != 'main.o')
command = ['gfortran', '-std=f2018', '-O2', '-ffree-line-length-none', '-fopenmp',
           '-I' + str(frozen / 'tmp/release'), '-J' + str(build), str(report / src.name),
           *map(str, objects), '-llapack', '-lblas', '-o', str(build / 'probe_interfaces')]
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
result = subprocess.run(command, capture_output=True, text=True, env=env)
(report / 'build.log').write_text(result.stdout + result.stderr)
assert result.returncode == 0, result.stderr
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
provenance = {
    'base_head': '9f2ed21668269bc60c8dbeaa5853b5700a081ac2',
    'operator_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=frozen, text=True).strip(),
    'operator_diff_from_base': subprocess.check_output(['git', 'diff', '8a85e37', '9f2ed21', '--', 'scr'], cwd=base, text=True),
    'command': command,
    'sources_sha256': {str(p.relative_to(frozen)): sha(p) for p in (frozen / 'scr').glob('*.f90')},
    'probe_source_sha256': sha(src), 'probe_binary_sha256': sha(build / 'probe_interfaces'),
    'operator_binary_sha256': sha(frozen / 'bin/prims_f'), 'runs': {},
}
assert not provenance['operator_diff_from_base']
for case in ('D_SS1_fixed', 'D_SS3_fixed'):
    target = report / case
    target.mkdir(exist_ok=True)
    original = base.parent / 'PRIMS-F-burr15-results-20261119/v1/runs' / case
    controls = json.loads((original / 'manifest.json').read_text())['controls']
    text = (original / 'run_metadata.txt').read_text()
    scales = text.split('Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u:')[1].splitlines()[0].split()
    values = [float(scales[-1]), controls['ps'], controls['jl'], controls['jg'], controls['mul'], controls['mug']]
    (target / 'probe_parameters.dat').write_text(' '.join(format(v, '.17e') for v in values) + ' 0\n')
    shutil.copyfile(original / 'last_rejected_state.dat', target / 'last_rejected_state.dat')
    result = subprocess.run([str(build / 'probe_interfaces')], cwd=target, env=env,
                            capture_output=True, text=True, timeout=120)
    (target / 'probe.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr
    provenance['runs'][case] = {'snapshot_sha256': sha(original / 'last_rejected_state.dat'),
                                'metadata_sha256': sha(original / 'run_metadata.txt'),
                                'operator_data_sha256': sha(target / 'state_operators.dat'),
                                'original': str(original), 'scales': list(map(float, scales)),
                                'controls': controls}
    print(case, result.stdout, flush=True)
(report / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
