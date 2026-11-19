from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119')
HERE = Path(__file__).resolve().parent
out = ROOT / 'active_set_smokes_v3_a179e51/SS3_terminal_80digits'
report_file = out / 'high_precision.json'
original = json.loads(report_file.read_text())
report = dict(original)
report['norm_semantics'] = 'Only physically valid assembled proposals have a residual norm. OLD was not assembled for this dt. Invalid assembly placeholders remain only in the raw CSV.'
for row in report['records']:
    if row['method'] == 'OLD' or row['physical_status_double'] != 0:
        row['norm_double'] = None
report['raw_report_sha256'] = hashlib.sha256(report_file.read_bytes()).hexdigest()
report['normalizer_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(out / 'high_precision_qualified.json').write_text(json.dumps(report, indent=2) + '\n')
for row in report['records']:
    if row['representation'] != 'exact_binary64':
        continue
    value = row['high_precision']
    print(row['method'], value['k'], 'Fend=' + value['endpoint_residual_m_s'],
          'bound=' + value['minimum_lower_bound_m_s'], row['classification'])

source = ROOT / 'active_set_smokes_v3_a179e51/D_SS3/output'
fixture = ROOT / 'offline_replay/D_SS3_active_terminal_v3'
fixture.mkdir(exist_ok=False)
shutil.copy2(source / 'last_rejected_state.dat', fixture / 'last_rejected_state.dat')
shutil.copy2(out / 'fixture/probe_parameters.dat', fixture / 'probe_parameters.dat')
# The independent frozen physical API is 69813ee. Its physical R is verified
# bitwise against this a179 snapshot before its double geometry is used.
from offline_probe import Probe, np
p = Probe('D_SS3_active_terminal_v3')
r = p.evaluate(p.trial)[0]
serialized = np.loadtxt(fixture / 'last_rejected_state.dat', skiprows=2)[:, 9:12].reshape(-1)
if not np.array_equal(r, serialized):
    raise RuntimeError('Frozen physical API residual differs from the a179 snapshot.')
completed = subprocess.run(['python3', str(HERE / 'inspect_terminal_geometry.py'),
                            '--case', 'D_SS3_active_terminal_v3'], text=True, capture_output=True)
(out / 'double_geometry.log').write_text(completed.stdout + completed.stderr)
completed.check_returncode()
for suffix in ('.json', '.npz'):
    shutil.copy2(ROOT / ('offline_replay/D_SS3_active_terminal_v3_linear_geometry' + suffix),
                 out / ('double_geometry' + suffix))
provenance = {'physical_api_build_manifest': str(ROOT / 'offline_replay/build_manifest.json'),
              'physical_api_build_manifest_sha256': hashlib.sha256((ROOT / 'offline_replay/build_manifest.json').read_bytes()).hexdigest(),
              'residual_bitwise_equal_to_a179_snapshot': True,
              'scope': 'Independent double geometry of the same unchanged physical residual; production direction replay uses a179 objects separately.'}
(out / 'double_geometry_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
