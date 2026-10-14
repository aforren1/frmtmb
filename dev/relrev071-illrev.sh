#!/usr/bin/env bash
# Release reviewer of 0.71.0: lane illse's reviewer's exact-reference
# scripts on rellib-r9 (the consolidation's copies in
# dev/rel071-log/reruns/r9/illse-rev/s-merged, which name rellib-r9),
# under ref and OpenBLAS 0.3.26 LAPACK (the consolidation ran 0.3.32),
# compared after the arm line with the consolidation's merged logs.
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
S=$REL/dev/rel071-log/reruns/r9/illse-rev
W=$REL/dev/relrev071-log/illse-rev
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$REL/dev/relrev071-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export NOT_CRAN=true REV_ARM=lane
rm -rf "$W"; mkdir -p "$W"
run() { # <cfg> <log> <script> args
  local c=$1 lg=$2 s=$3; shift 3
  local R="$RS"
  if [ $c = ob26 ]; then R="$OB"; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
  else unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS; fi
  "$R" "$S/s-merged/$s" "$@" < /dev/null > "$W/$lg-$c.txt" 2>&1 &
}
cd "$REL"
for c in ref ob26; do
  run $c lmm8 illse-rev-lmm.R 101:103 8
  run $c glm-binom illse-rev-glm.R binom 201:204
  run $c edge4 illse-rev-edge4.R
done
run ref lmm7 illse-rev-lmm.R 101:103 7
run ref insight illse-rev-insight.R
run ref edge8 illse-rev2-edge8.R
wait
for f in "$W"/*-ref.txt; do
  b=$(basename "$f" -ref.txt)
  g="$S/$b-ref-merged.txt"
  if diff -q <(tail -n +2 "$f") <(tail -n +2 "$g") > /dev/null; then
    echo "SAME as the consolidation's ref log: $b"
  else
    echo "DIFFERENT from the consolidation's ref log: $b ($(diff <(tail -n +2 "$f") <(tail -n +2 "$g") | grep -c '^[<>]') lines)"
  fi
done > "$W/summary.txt"
