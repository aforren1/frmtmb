#!/usr/bin/env bash
# Fix round: a job list of "<package> <test file>" lines over the BLAS
# settings the coordinator named, at most 3 R processes at once.
#   bash dev/ci082-matrix2.sh <lib|base> <outdir> <job list>
# Settings: reference BLAS; OpenBLAS 0.3.26 and 0.3.32 with LAPACK, with
# OPENBLAS_CORETYPE unset, Haswell, Sandybridge. OPENBLAS_NUM_THREADS=4 as
# on the runners, except test-se-limit.R, which runs 1 thread: under a
# full machine 4 threads stalled its nonlinear gp() fit for over 15
# minutes, base and lane alike (dev/ci082-findings.md).
set -u
ARM="$1"; OUT="$2"; LIST="$3"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export NOT_CRAN=true
mkdir -p "$OUT"
n=0
while read -r P F; do
  [ -n "$P" ] || continue
  for CFG in ref 0.3.26-default 0.3.26-Haswell 0.3.26-Sandybridge \
             0.3.32-default 0.3.32-Haswell 0.3.32-Sandybridge; do
    if [ "$CFG" = ref ]; then
      R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"; CT=default
    else
      R="dev/ciharden-out/Rob${CFG%%-*}-lapack/bin/Rscript.exe"; CT="${CFG#*-}"
    fi
    TH=4
    [ "$(basename "$F")" = test-se-limit.R ] && TH=1
    LOG="$OUT/$P--$(basename "$F" .R)--$CFG.txt"
    (
      export OPENBLAS_NUM_THREADS=$TH
      if [ "$CT" = default ]; then unset OPENBLAS_CORETYPE; else
        export OPENBLAS_CORETYPE="$CT"; fi
      "$R" dev/ci082-run1.R "$ARM" "$P" "$F" > "$LOG" 2>&1
      echo "exit $?" >> "$LOG"
    ) &
    n=$((n + 1))
    if [ "$n" -ge "${MAXJ:-3}" ]; then wait -n; n=$((n - 1)); fi
  done
done < "$LIST"
wait
grep -a -H "^RESULT\|^exit [1-9]" "$OUT"/*.txt
echo "MATRIX DONE"
