#!/bin/bash
# P0.4 isolated benchmark: ONE process at a time, pinned to one core, 1 thread.
R=$(cd "$(dirname "$0")/.." && pwd); B=$R/bench
S5=~/Documentos/PRIMS-F-chexal-riser-results-20261121/bin/prims_f_seal5_46e8830
J1=$(ls $R/bin/prims_f_jac1_* | head -1)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1
{ echo "date $(date -Is)"; uname -a; lscpu | grep -E "Model name|MHz|Thread|Core"; cat /sys/devices/system/cpu/cpu3/cpufreq/scaling_governor 2>/dev/null;
  gfortran --version | head -1; ldd $J1 | grep -iE "blas|lapack|mkl"; sha256sum $S5 $J1; } > $B/environment.txt 2>&1
for run in ss1_s5 ss1_A ss1_C16fd ss1_C16ij rep_s5 rep_C16fd rep_C16ij; do
  bin=$J1; case $run in *_s5) bin=$S5;; esac
  rm -rf $B/runs/$run; mkdir -p $B/runs/$run
  while pgrep -f "bin/prims_f" >/dev/null; do sleep 30; done
  ( cd $B/runs && /usr/bin/time -v taskset -c 3 timeout 10800 $bin $B/decks/$run.nml > $run.stdout 2> $run.time )
  echo "$run rc=$? $(grep -E 'termination_reason|reached_time' $B/runs/$run/run_status.txt | tr '\n' ' ')"
done
echo bench done
