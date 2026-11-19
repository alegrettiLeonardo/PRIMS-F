import ctypes as ct
import json
import subprocess
from pathlib import Path
from offline_probe import Probe,ROOT,np,ND,IP
here=Path(__file__).resolve().parent
build=ROOT/'build'
command=['gfortran','-O2','-fPIC','-shared','-fopenmp','-I'+str(build),str(here/'cell_api.f90'),
         str(build/'offline_api.so'),'-Wl,-rpath,'+str(build),'-o',str(build/'cell_api.so')]
subprocess.run(command,check=True,capture_output=True)
p=Probe('D_SS3_active_terminal')
lib=ct.CDLL(str(build/'cell_api.so'))
lib.probe_cell.argtypes=[ND,ct.c_int,ND,IP]
records={}
for label,x in [('accepted_old',p.old),('linearization',p.trial)]:
    data=np.empty(15);status=ct.c_int()
    lib.probe_cell(x,69,data,ct.byref(status))
    records[label]={'status':status.value,'values':data.tolist()}
out={'columns':['alpha','rhol_kg_m3','rhog_kg_m3','p_Pa','ul_m_s','ug_m_s','ml_kg_m3','mg_kg_m3','M_kg_m2_s',
                'theta_from_horizontal_rad','D_m','g_m_s2','mu_l_Pa_s','mu_g_Pa_s','sigma_N_m'],'records':records}
(ROOT/'D_SS3_active_terminal_cell.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
