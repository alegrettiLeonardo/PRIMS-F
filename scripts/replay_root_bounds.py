"""Reconstruct exact rejected states and check SS3 against an independent bound."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
from generate_bhagwat_fixtures import coefficients


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def freeze_probe(output, objects):
    source = (ROOT/"validation/diagnostics/burr15_rejection_replay_20261119/probe_rejection_v2.f90").read_text()
    assert source.count("    if(profile==2.and.code==0)write") == 1
    source = source.replace("    if(profile==2.and.code==0)write", """    if(profile==2.and.code==0)then
      open(newunit=iu,file='linearization.csv',status='new')
      write(iu,'(a)')'cell,r1,r2,r3,A11,A21,A31,A12,A22,A32,A13,A23,A33,'// &
        'B11,B21,B31,B12,B22,B32,B13,B23,B33,C11,C21,C31,C12,C22,C32,C13,C23,C33'
      do i=3,104
        write(iu,'(i0,30(a,es26.17e3))')i,(',',r(k,i-2),k=1,3), &
          ((',',ab(k,j,i-2),k=1,3),j=1,3),((',',bb(k,j,i-2),k=1,3),j=1,3), &
          ((',',cb(k,j,i-2),k=1,3),j=1,3)
      enddo
      close(iu)
      open(newunit=iu,file='proposal_cells.csv',status='new')
      write(iu,'(a)')'k,cell,lambda,theta,U1,U2,U3'
      do i=53,69
        if(i/=53.and.i/=54.and.i/=69)cycle
        write(iu,'(i0,a,i0,5(a,es26.17e3))')-1,',',i,',',0.0_wp,',',cfg%thetav(i), &
          (',',state(i,k),k=1,3)
      enddo
      close(iu)
    endif
    if(profile==2.and.code==0)write""")
    # The existing program declares k, i and state. j is used only to serialize blocks.
    source = source.replace("  implicit none", "  implicit none\n  integer :: j")
    needle = "        call sim%assemble(candidate,old,dt,r,ab,bb,cb,nt,fin,fout,code,.false.)"
    assert source.count(needle) == 1
    prefix = """        open(newunit=iu,file='proposal_cells.csv',status='old',position='append')
        do i=53,69
          if(i/=53.and.i/=54.and.i/=69)cycle
          write(iu,'(i0,a,i0,5(a,es26.17e3))')k,',',i,',',lambda,',',cfg%thetav(i), &
            ',',candidate(i,1),',',candidate(i,2),',',candidate(i,3)
        enddo
        close(iu)
"""
    source = source.replace(needle, prefix+needle)
    probe = output/"probe_limits.f90"
    probe.write_text(source)
    binary = output/"probe_limits"
    selected = sorted(p for p in objects.glob("*.o") if p.name != "main.o")
    if not selected:
        raise ValueError("No compiled Fortran objects")
    command = ["gfortran","-std=f2018","-O2","-ffree-line-length-none","-fopenmp",
               "-I"+str(objects),str(probe),*map(str,selected),"-llapack","-lblas","-o",str(binary)]
    subprocess.run(command, check=True)
    return binary, {"command":command,"probe_sha256":digest(probe),"binary_sha256":digest(binary),
                    "objects_sha256":{obj.name:digest(obj) for obj in selected}}


def bound(row, controls, velocity_scale):
    """SI residual F(q) = C0*j + Ud - ug, q = gas momentum fraction."""
    m = mp.mpf
    ml, mg = m(row["U1"])*m("997.05"), m(row["U2"])*m("997.05")
    total = m(row["U3"])*m("997.05")*m(velocity_scale)
    angle = m(row["theta"])
    cl, cg, p0, rho0 = (m(str(controls[k])) for k in ("cl","cg_fixed","p_l0","rho_l0"))
    a = rho0-p0/cl**2
    q = mg*cg**2/cl**2
    b = a-ml-q
    pressure = (mp.sqrt(b*b+4*a*q)-b)*cl**2/2
    rg, rl = pressure/cg**2, rho0+(pressure-p0)/cl**2
    alpha = mg/rg
    d, eps, mu, sigma = (m(str(controls[k])) for k in ("d","eps","mu_l","sigma"))
    ugmax, ulmax = total/mg, total/ml
    j0, j1 = (1-alpha)*ulmax, alpha*ugmax

    def residual(qgas):
        ul, ug = (1-qgas)*ulmax,qgas*ugmax
        c0, drift, _ = coefficients((alpha,rl,rg,ul,ug,angle,d,eps,mu,sigma))
        return c0*((1-qgas)*j0+qgas*j1)+drift-ug

    c0, drift, _ = coefficients((alpha,rl,rg,0,ugmax,angle,d,eps,mu,sigma))
    cut = min(m(".5"),drift/(2*ugmax))
    rough = eps/(m("3.7")*d)
    min_re = rl*j0*d/mu
    cut_re = rl*((1-cut)*j0+cut*j1)*d/mu
    hthird = m(1)/3+4*mp.log10(rough+m("1.256")/(3*min_re))
    hone = 1+4*mp.log10(rough+m("1.256")/cut_re)
    minimum_bound = min(drift-cut*ugmax, drift+min(m(0),alpha*c0-1)*ugmax)
    end = residual(m(1))
    conditions = (0<=angle<=mp.pi/2 and 0<alpha<1 and rl>rg>0 and j1>=j0>0
                  and cut>0 and hthird<0 and hone<0)
    # Derivative from inside the admissible interval, including x=1 at its endpoint.
    step = m("1e-24")
    endpoint_slope = (end-residual(1-step))/step
    values = {"pressure_Pa":pressure, "alpha":alpha,"rho_l_kg_m3":rl,"rho_g_kg_m3":rg,
              "u_g_max_m_s":ugmax,"minimum_lower_bound_m_s":minimum_bound,
              "endpoint_residual_m_s":end,"endpoint_dF_dq_m_s":endpoint_slope,
              "split_gas_momentum_fraction":cut,"H_at_one_third_low":hthird,"H_at_one_high":hone}
    root = None
    if end<0 and residual(cut)>0:
        lo,hi = cut,m(1)
        for _ in range(115):
            mid = (lo+hi)/2
            if residual(mid)>0:
                lo=mid
            else:
                hi=mid
        root = (lo+hi)/2
        values.update(root_q=root,root_ul_m_s=(1-root)*ulmax,root_ug_m_s=root*ugmax)
    return {"k":int(row["k"]), "lambda":row["lambda"].strip(),
            "excluded":bool(conditions and minimum_bound>0),
            **{name:mp.nstr(value,35) for name,value in values.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--objects",type=Path,required=True)
    parser.add_argument("--runs",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True)
    binary, report = freeze_probe(out,args.objects.resolve())
    report["git_head"] = subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
    report["script_sha256"] = digest(Path(__file__))
    env = dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    report["cases"] = {}
    for name in ("D_SS1_fixed","D_SS3_fixed"):
        previous = args.runs/name
        manifest = json.loads((previous/"manifest.json").read_text())
        controls = manifest["controls"]
        scales = (previous/"run_metadata.txt").read_text().split(
            "Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u:")[1].splitlines()[0].split()
        case = {"snapshot_sha256":digest(previous/"last_rejected_state.dat"),"modes":{}}
        for flag in (0,1):
            folder = out/f"{name}_{flag}"
            folder.mkdir()
            shutil.copy2(previous/"last_rejected_state.dat",folder/"last_rejected_state.dat")
            values = [scales[-1],*(controls[k] for k in ("ps","jl","jg","mul","mug")),flag]
            (folder/"probe_parameters.dat").write_text(" ".join(map(str,values))+"\n")
            captured = subprocess.run([binary],cwd=folder,env=env,capture_output=True,text=True,check=True)
            (folder/"probe.log").write_text(captured.stdout)
            rows = []
            if name == "D_SS3_fixed":
                with (folder/"proposal_cells.csv").open() as stream:
                    rows = [bound(row,controls,scales[-1]) for row in csv.DictReader(stream) if row["cell"]=="69"]
            case["modes"][flag] = {"bounds_cell69":rows,
                                   "log_sha256":digest(folder/"probe.log"),
                                   "proposal_cells_sha256":digest(folder/"proposal_cells.csv"),
                                   "linearization_sha256":digest(folder/"linearization.csv")}
        report["cases"][name] = case
    (out/"interval_replay.json").write_text(json.dumps(report,indent=2)+"\n")


if __name__ == "__main__":
    main()
