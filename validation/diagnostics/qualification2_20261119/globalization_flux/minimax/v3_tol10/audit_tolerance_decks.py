from pathlib import Path
import hashlib
import json
import re

base = Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119')
reference = base / 'active_set_smokes_v3_a179e51'
target = base / 'active_set_smokes_v3_tol10_a179e51'
report = {'scope': 'Identical frozen operator, only tol_res tightened; out_dir differs by cohort.', 'cases': []}
for case in ('D_SS1', 'D_SS3'):
    old = reference / case / (case + '.nml')
    new = target / case / (case + '.nml')
    before = re.sub(r"out_dir\s*=\s*'[^']*'", "out_dir='COHORT'", old.read_text())
    after = re.sub(r"out_dir\s*=\s*'[^']*'", "out_dir='COHORT'", new.read_text())
    expected, count = re.subn(r'\btol_res\s*=\s*[0-9.eEdD+-]+', 'tol_res=1e-10', before)
    if count != 1 or expected != after:
        raise RuntimeError('Unexpected difference beyond tol_res/out_dir in ' + case)
    if not old.read_bytes().endswith(b'\n') or not new.read_bytes().endswith(b'\n'):
        raise RuntimeError('Missing final namelist newline.')
    report['cases'].append({'case': case, 'verified': True,
        'reference_deck_sha256': hashlib.sha256(old.read_bytes()).hexdigest(),
        'tight_deck_sha256': hashlib.sha256(new.read_bytes()).hexdigest(),
        'tol_res_reference': '1e-08', 'tol_res_tight': '1e-10',
        'tol_inc_unchanged': re.search(r'\btol_inc\s*=\s*([0-9.eEdD+-]+)', after).group(1)})
(target / 'deck_audit.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
