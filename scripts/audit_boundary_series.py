#!/usr/bin/env python3
"""Replay inlet/outlet characteristics from stored accepted field frames."""
import argparse
import csv
import hashlib
from itertools import groupby
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from validation_evidence import metadata


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()


def write_input(run,path):
    p=metadata(run/"run_metadata.txt")
    # Values follow the exact scaling transformation used by the driver.
    line=(run/"run_metadata.txt").read_text().split("Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u:")[1].splitlines()[0]
    _,_,_,_,_,wp,wr,wu=map(float,line.split())
    area=np.pi*p["d"]**2/4
    values=[wu,wr,wp,p["d"],area,p["eps"],p["g"],p["mu_l"],p["mu_g"],p["sigma"],
            p["cl"]*np.sqrt(wr/wp),p["cg_fixed"]*np.sqrt(wr/wp),p["rho_l0"]/wr,
            p["p_l0"]/wp,p["ps"]/wp,p["mul"]/(wr*wu),p["mug"]/(wr*wu),p["jl"]/wu]
    count=0
    with path.open("x") as out,(run/"fields.csv").open() as stream:
        out.write(" ".join(format(v,".17e") for v in values)+"\n")
        for t,group in groupby(csv.DictReader(stream),key=lambda r:float(r["time_s"])):
            first=next(group)
            last=first
            for row in group:
                last=row
            for side,row in enumerate((first,last),1):
                state=[t,side,float(row["theta_rad"]),float(row["alpha"]),float(row["rhol_kg_m3"])/wr,
                       float(row["rhog_kg_m3"])/wr,float(row["p_Pa"])/wp,
                       float(row["ul_m_s"])/wu,float(row["ug_m_s"])/wu]
                out.write(f"{t:.17e} {side} "+" ".join(format(v,".17e") for v in state[2:])+"\n")
            count+=1
    return count


def audit(path,expected_frames):
    result={}
    with path.open() as stream:
        for row in csv.DictReader(stream):
            key=f"{int(row['side'])}_{int(row['kind'])}"
            st=int(row["status"])
            bucket=result.setdefault(key,{"samples":0,"status_counts":{},"incoming_counts":{},
                                         "eigen_min_m_s":np.inf,"eigen_max_m_s":-np.inf,
                                         "max_imaginary_m_s":0.,"min_constraint_singular_value":np.inf,
                                         "negative_ul_samples":0,"negative_ug_samples":0})
            bucket["samples"]+=1
            bucket["status_counts"][st]=bucket["status_counts"].get(st,0)+1
            if st:
                continue
            eigen=np.array([float(row[f"e{k}_m_s"]) for k in range(1,4)])
            imag=np.array([float(row[f"i{k}_m_s"]) for k in range(1,4)])
            tolerance=64*np.finfo(float).eps*max(1.,np.max(np.abs(eigen)))
            incoming=eigen>tolerance
            if int(row["side"])==2:
                incoming=eigen< -tolerance
            n=int(np.count_nonzero(incoming))
            bucket["incoming_counts"][n]=bucket["incoming_counts"].get(n,0)+1
            bucket["eigen_min_m_s"]=min(bucket["eigen_min_m_s"],float(min(eigen)))
            bucket["eigen_max_m_s"]=max(bucket["eigen_max_m_s"],float(max(eigen)))
            bucket["max_imaginary_m_s"]=max(bucket["max_imaginary_m_s"],float(max(abs(imag))))
            bucket["negative_ul_samples"]+=float(row["ul_m_s"])<0
            bucket["negative_ug_samples"]+=float(row["ug_m_s"])<0
            projection=np.array([[float(row[f"c{i}{j}"]) for j in range(1,4)] for i in range(1,4)])
            constraints=projection[:2,incoming] if int(row["side"])==1 else projection[2:,incoming]
            if constraints.shape[0]==n and n>0:
                smallest=float(np.linalg.svd(constraints,compute_uv=False).min())
                bucket["min_constraint_singular_value"]=min(bucket["min_constraint_singular_value"],smallest)
    if set(result)!={"1_1","1_2","2_1","2_2"} or any(bucket["samples"]!=expected_frames for bucket in result.values()):
        raise ValueError("Incomplete boundary characteristic coverage")
    for bucket in result.values():
        for key,value in bucket.items():
            if isinstance(value,float) and not np.isfinite(value):
                bucket[key]=None
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run",type=Path,required=True)
    parser.add_argument("--objects",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]
    source=root/"validation/diagnostics/qualification2_20261119/probe_boundary_series.f90"
    objects=sorted(p for p in args.objects.glob("*.o") if p.name!="main.o")
    command=["gfortran","-std=f2018","-O2","-ffree-line-length-none","-fopenmp","-I"+str(args.objects),
             str(source),*map(str,objects),"-llapack","-lblas","-o",str(out/"probe")]
    subprocess.run(command,check=True,capture_output=True)
    count=write_input(args.run,out/"boundary_input.dat")
    env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    subprocess.run([out/"probe",out/"boundary_input.dat",out/"boundary_characteristics.csv"],check=True,env=env)
    report={"run":str(args.run),"input_frames":count,"boundary_summary":audit(out/"boundary_characteristics.csv",count),
            "provenance":{str(path):sha(path) for path in (args.run/"fields.csv",args.run/"manifest.json",source,Path(__file__))},
            "command":command,
            "scope":"Accepted field sampling only; rows reconstruct conservative states from printed primitives. "
                    "Interior and imposed states share the documented branch. Unsampled Newton proposals are separate."}
    (out/"summary.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
