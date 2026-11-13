#!/usr/bin/env bash
# The lane's four test files over every BLAS setting tried, at most 3 R
# processes at once.
#   bash dev/ci082-matrix.sh <lib|base> <outdir> <file...>
# Settings: reference BLAS; OpenBLAS 0.3.26 and 0.3.32 (LAPACK too) with
# OPENBLAS_CORETYPE unset, Haswell, Zen, Sandybridge (SkylakeX needs
# AVX-512, which this CPU lacks). OPENBLAS_NUM_THREADS=4 as on the runners,
# or OB_THREADS: under a full machine, 4 threads stalled the nonlinear gp()
# fit of test-se-limit.R for over 15 minutes, base and lane alike.
set -u
ARM="$1"; OUT="$2"; shift 2
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export NOT_CRAN=true OPENBLAS_NUM_THREADS="${OB_THREADS:-4}"
mkdir -p "$OUT"
n=0
for F in "$@"; do
  for CFG in ref 0.3.26-default 0.3.26-Haswell 0.3.26-Zen 0.3.26-Sandybridge \
             0.3.32-default 0.3.32-Haswell 0.3.32-Zen 0.3.32-Sandybridge; do
    if [ "$CFG" = ref ]; then
      R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"; CT=default
    else
      R="dev/ciharden-out/Rob${CFG%%-*}-lapack/bin/Rscript.exe"; CT="${CFG#*-}"
    fi
    LOG="$OUT/$(basename "$F" .R)-$CFG.txt"
    (
      if [ "$CT" = default ]; then unset OPENBLAS_CORETYPE; else
        export OPENBLAS_CORETYPE="$CT"; fi
      "$R" dev/ci082-run1.R "$ARM" frmtmb "$F" > "$LOG" 2>&1
      echo "exit $?" >> "$LOG"
    ) &
    n=$((n + 1))
    if [ "$n" -ge 3 ]; then wait -n; n=$((n - 1)); fi
  done
done
wait
grep -H "RESULT\|^exit" "$OUT"/*.txt
echo "MATRIX DONE"
