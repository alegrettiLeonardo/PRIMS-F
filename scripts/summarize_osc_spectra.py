#!/usr/bin/env python3
"""One common finite-window spectral definition for nominal and BC sensitivities."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
from plot_osc_experimental import digitized_pressure,spectrum
from validation_evidence import load_runs,time_statistics


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()


def spectral_peaks(time,signal):
    time,signal=np.asarray(time),np.asarray(signal)
    if len(time)<2 or len(signal)!=len(time) or not np.isfinite(signal).all() or not np.isfinite(time).all():
        raise ValueError('Finite pressure samples and at least two timestamps are required')
    if np.any(np.diff(time)<=0):
        raise ValueError('Time must increase')
    excursion=float(np.ptp(signal))
    useful=excursion>=100.
    if excursion==0:
        return {"frequency_resolution_Hz":float((len(time)-1)/(len(time)*(time[-1]-time[0]))),
                "three_largest_local_peaks":[],"ratio_two_strongest_frequency_bins":None,
                "window_excursion_Pa":0.,"window_excursion_at_least_100Pa":False,
                "peaks_interpretation":"NO_SPECTRAL_VARIATION; no period inferred"}
    frequency,power=spectrum(time,signal)
    peaks=[index for index in find_peaks(power)[0] if 0<frequency[index]<=.1]
    selected=sorted(peaks,key=lambda index:power[index],reverse=True)[:3]
    step=float(frequency[1]-frequency[0])
    records=[{"frequency_bin_Hz":float(frequency[index]),"relative_power":float(power[index]),
              "bin_half_width_Hz":step/2} for index in selected]
    ratio=None
    if len(records)>=2:
        ordered=sorted(records[:2],key=lambda row:row["frequency_bin_Hz"])
        ratio=ordered[1]["frequency_bin_Hz"]/ordered[0]["frequency_bin_Hz"]
    return {"frequency_resolution_Hz":step,"three_largest_local_peaks":records,
            "ratio_two_strongest_frequency_bins":ratio,
            "window_excursion_Pa":excursion,"window_excursion_at_least_100Pa":useful,
            "peaks_interpretation":("Finite-window frequency bins; cycle qualification is separate" if useful else
                                    "NO_USEFUL_WINDOW_VARIATION; normalized peaks do not qualify a period"),
            "window":"Hann, mean removed, relative power normalized separately for each series"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--digitization",type=Path,required=True)
    parser.add_argument("runs",type=Path,nargs="+")
    args=parser.parse_args()
    te,pe,lo,hi=digitized_pressure(args.digitization)
    result={"experiment":spectral_peaks(te,pe),
            "digitization_sha256":sha(args.digitization),
            "interpretation":"Bins of a finite published excerpt. A ratio near two is compatible with a harmonic "
                             "or subharmonic structure; the largest peak alone does not identify a physical mode. "
                             "Calibration, curve thickness and nonstationarity add uncertainty beyond bin width.",
            "runs":[]}
    duration=te[-1]-te[0]
    for folder in args.runs:
        record,=load_runs(folder.parent,pattern=folder.name)
        data=record["data"]
        stop=record["reached_time_s"]
        if not record["completed"] or stop<duration:
            result["runs"].append({"run":record["run"],"status":"INCOMPLETE"})
            continue
        start=stop-duration
        time=np.array(data["time_s"])
        mask=time>=start
        signal=np.array(data["p_junction_Pa"],dtype=float)
        values=spectral_peaks(time[mask],signal[mask])
        values.update(run=record["run"],status="COMPLETED",start_s=float(time[mask][0]),end_s=float(stop),
                      pressure_column="p_junction_Pa",statistics=time_statistics(data,"p_junction_Pa",start,stop),
                      manifest_sha256=sha(folder/"manifest.json"),monitor_sha256=sha(folder/"monitor.csv"))
        result["runs"].append(values)
    result["script_sha256"]=sha(Path(__file__))
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")


if __name__=="__main__":
    main()
