from pathlib import Path
import subprocess
import json
import sys
HERE=Path(__file__).resolve().parent
OUT=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119/offline_replay')
for case in ('D_SS1','D_SS3'):
    for factor in (0.25,0.5,2,4):
        result=subprocess.run([sys.executable,str(HERE/'active_replay.py'),'--case',case,'--fd-factor',str(factor)],text=True,capture_output=True)
        (OUT/(case+f'_fd{factor:g}.log')).write_text(result.stdout+result.stderr)
        report=json.loads((OUT/(case+f'_active_report_fd{factor:g}.json')).read_text())
        print(case,factor,result.returncode,report['termination'],report['final_norm'],report['final_margin_si'],flush=True)
