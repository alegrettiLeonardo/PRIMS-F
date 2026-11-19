from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

repo = Path('/home/leonardoab/Documentos/PRIMS-F-numerics-20261119')
frozen = repo.parent / 'PRIMS-F-burr15-regression-8a85e37-20261119'
report = repo / 'validation/diagnostics/numerics_p0p1_20261119'
target = report / 'SS1_minimal_replay'
target.mkdir(exist_ok=True)
build = repo / 'tmp/p0p1_probe'
source = repo / 'validation/diagnostics/burr15_rejection_replay_20261119/probe_rejection_v2.f90'
objects = sorted(p for p in (frozen / 'tmp/release').glob('*.o') if p.name != 'main.o')
command = ['gfortran', '-std=f2018', '-O2', '-ffree-line-length-none', '-fopenmp',
           '-I' + str(frozen / 'tmp/release'), '-J' + str(build), str(source),
           *map(str, objects), '-llapack', '-lblas', '-o', str(build / 'replay_ss1')]
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
process = subprocess.run(command, capture_output=True, text=True, env=env)
assert process.returncode == 0, process.stderr
shutil.copyfile(report / 'D_SS1_fixed/last_rejected_state.dat', target / 'last_rejected_state.dat')
parameters = (report / 'D_SS1_fixed/probe_parameters.dat').read_text().split()
parameters[-1] = '1'
(target / 'probe_parameters.dat').write_text(' '.join(parameters) + '\n')
process = subprocess.run([str(build / 'replay_ss1')], cwd=target, env=env,
                         capture_output=True, text=True, timeout=120)
(target / 'replay.log').write_text(process.stdout + process.stderr)
assert process.returncode == 0, process.stderr
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
record = {'base_head': '9f2ed21668269bc60c8dbeaa5853b5700a081ac2',
          'operator_head': '8a85e37319995b37795d2c798d6f865b863b1834',
          'source_unchanged_between_heads': True, 'build_command': command,
          'returncode': process.returncode, 'source_sha256': sha(source),
          'binary_sha256': sha(build / 'replay_ss1'),
          'snapshot_sha256': sha(target / 'last_rejected_state.dat'),
          'log_sha256': sha(target / 'replay.log')}
lines = process.stdout.splitlines()
record['residual_max_abs_difference'] = float(next(line.split()[-1] for line in lines if line.startswith('residual difference')))
record['proposal_count'] = sum(line.startswith('proposal k,') for line in lines)
record['first_face_failures'] = sorted(set(line for line in lines if line.startswith('first face failure')))
record['proposal_statuses'] = sorted(set(int(line.split()[-3]) for line in lines if line.startswith('proposal k,')))
assert record['residual_max_abs_difference'] == 0
assert record['proposal_count'] == 18 and record['proposal_statuses'] == [3012]
(target / 'provenance.json').write_text(json.dumps(record, indent=2) + '\n')
shutil.copyfile(Path(__file__), report / Path(__file__).name)
print(json.dumps({key: value for key, value in record.items() if key != 'build_command'}, indent=2))
