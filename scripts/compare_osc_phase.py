#!/usr/bin/env python3
"""Compare settled accepted OSC cycles, means, phase waveforms and integral balances."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from analyze_osc_envelopes import pressure_envelope,last_ten_cycles
from validation_evidence import load_runs,time_statistics

CHANNELS=("p_riser_base_Pa","p_junction_Pa","holdup_liquid_base","holdup_liquid_out",
          "mdot_l_base_kg_s","mdot_g_base_kg_s","mdot_l_out_kg_s","mdot_g_out_kg_s",
          "Ml_kg","Mg_kg","Ml_pipe_kg","Mg_pipe_kg","Ml_riser_kg","Mg_riser_kg")


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()


def phase_waveform(time,signal,bounds,points=256):
    phase=np.arange(points)/points
    cycles=np.array([np.interp(left+(right-left)*phase,time,signal)
                     for left,right in zip(bounds,bounds[1:])])
    mean=cycles.mean(axis=0)
    excursion=np.ptp(mean)
    marker=None
    if excursion>0:
        idx=int(np.argmax(mean))
        marker=float((phase[idx]+.5)%1-.5)
    return {"mean":mean.tolist(),"phase_at_maximum_relative_to_pressure_peak":marker,
            "phase_cycle_scatter_rms":float(np.sqrt(np.mean((cycles-mean)**2))),
            "excursion_of_phase_mean":float(excursion)}


def summarize(folder):
    run,=load_runs(folder.parent,pattern=folder.name)
    out={k:run[k] for k in ("run","completed","reached_time_s","provenance")}
    manifest=json.loads((folder/"manifest.json").read_text())
    out["controls"]=manifest["controls"]
    out["freeze"]=manifest["git_head"]
    out["binary_sha256"]=manifest["binary_sha256"]
    data=run["data"]
    if not run["completed"]:
        out["status"]="INCOMPLETE"
        return out
    t=np.array(data["time_s"])
    p=np.array(data["p_junction_Pa"])
    metrics=last_ten_cycles(t,p,pressure_envelope(t,p))
    out["cycles"]=metrics
    out["pressure_marker_column"]="p_junction_Pa"
    if not metrics["available"] or not metrics["useful_above_declared_detection_thresholds"]:
        out["status"]="NO_CYCLES"
        return out
    bounds=[entry["peak_time_s"] for entry in metrics["peak_records"]]
    out["window_start_s"]=bounds[0];out["window_end_s"]=bounds[-1]
    out["statistics"]={}
    out["waveforms"]={}
    drift=[]
    for key in CHANNELS:
        signal=data[key]
        if any(v is None for v in signal):
            raise ValueError(f"Missing {key}")
        stats=time_statistics(data,key,bounds[0],bounds[-1])
        out["statistics"][key]=stats
        out["waveforms"][key]=phase_waveform(t,signal,bounds)
        cycles=[time_statistics(data,key,a,b) for a,b in zip(bounds,bounds[1:])]
        for variable in ("mean","peak_to_peak"):
            first=np.mean([cycle[variable] for cycle in cycles[:5]])
            last=np.mean([cycle[variable] for cycle in cycles[5:]])
            scale=abs(first)
            if scale>0:
                drift.append(abs(float(last-first)/scale))
    out["largest_cycle_drift"]=max(drift)
    out["status"]="SETTLED" if metrics["amplitude_change_within_one_percent"] and max(drift)<=.01 else "UNSETTLED"
    out["hashes"]={name:sha(folder/name) for name in ("manifest.json","monitor.csv","steps.csv","fields.csv")}
    return out


def compare(left,right):
    if left["freeze"]!=right["freeze"] or left["binary_sha256"]!=right["binary_sha256"]:
        raise ValueError("Different operator versions")
    changes={key:{"left":left["controls"].get(key),"right":right["controls"].get(key)}
             for key in left["controls"].keys()|right["controls"].keys()
             if key!="out_dir" and left["controls"].get(key)!=right["controls"].get(key)}
    out={"left":left["run"],"right":right["run"],"both_settled":left["status"]==right["status"]=="SETTLED",
         "control_differences":changes}
    if "waveforms" not in left or "waveforms" not in right:
        return out
    for name in ("mean_period_s","mean_semi_amplitude_Pa"):
        first,last=left["cycles"][name],right["cycles"][name]
        relative=(last-first)/abs(first)
        out[name]={"left":first,"right":last,"relative_change":relative,"passes_one_percent":abs(relative)<=.01}
    out["channels"]={}
    for key in CHANNELS:
        a=np.array(left["waveforms"][key]["mean"])
        b=np.array(right["waveforms"][key]["mean"])
        excursion=float(np.ptp(a))
        error=b-a
        waveform={"rms_absolute_difference":float(np.sqrt(np.mean(error**2))),
                  "maximum_absolute_difference":float(np.max(abs(error))),
                  "rms_relative_to_left_excursion":float(np.sqrt(np.mean(error**2))/excursion) if excursion>0 else None,
                  "mean_difference":right["statistics"][key]["mean"]-left["statistics"][key]["mean"],
                  "maximum_phase_difference_cycles":None}
        am=left["waveforms"][key]["phase_at_maximum_relative_to_pressure_peak"]
        bm=right["waveforms"][key]["phase_at_maximum_relative_to_pressure_peak"]
        if am is not None and bm is not None:
            waveform["maximum_phase_difference_cycles"]=float((bm-am+.5)%1-.5)
        out["channels"][key]=waveform
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("runs",type=Path,nargs="+")
    args=parser.parse_args()
    result={"schema":"PRIMS-F.phase_comparison.v1","phase_grid_points":256,
            "phase_origin":"enclosed pressure peaks; no waveform shifting optimized per channel",
            "runs":[summarize(folder) for folder in args.runs]}
    reference=result["runs"][0]
    result["comparisons"]=[compare(reference,run) for run in result["runs"][1:]]
    result["script_sha256"]=sha(Path(__file__))
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")


if __name__=="__main__":
    main()
