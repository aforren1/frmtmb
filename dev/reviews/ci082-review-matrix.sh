#!/usr/bin/env bash
# Review of lane ci082: test files on the lane library over BLAS
# settings, at most 2 R processes at once (the third slot is kept for
# probes).
#   bash dev/reviews/ci082-review-matrix.sh <outdir> <cfgs> <file...>
# cfgs: comma list of ref, 0.3.26-default, 0.3.26-Haswell, ...
set -u
OUT="$1"; CFGS="$2"; shift 2
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export NOT_CRAN=true OPENBLAS_NUM_THREADS="${OB_THREADS:-4}"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
mkdir -p "$OUT"
n=0
for F in "$@"; do
  for CFG in ${CFGS//,/ }; do
    if [ "$CFG" = ref ]; then
      R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"; CT=default
    else
      R="dev/ciharden-out/Rob${CFG%%-*}-lapack/bin/Rscript.exe"; CT="${CFG#*-}"
    fi
    LOG="$OUT/$(basename "$F" .R)-$CFG.txt"
    (
      if [ "$CT" = default ]; then unset OPENBLAS_CORETYPE; else
        export OPENBLAS_CORETYPE="$CT"; fi
      "$R" dev/ci082-run1.R lib frmtmb "$F" > "$LOG" 2>&1
      echo "exit $?" >> "$LOG"
    ) &
    n=$((n + 1))
    if [ "$n" -ge "${JOBS:-2}" ]; then wait -n; n=$((n - 1)); fi
  done
done
wait
grep -H "RESULT\|^exit\|^lib:" "$OUT"/*.txt
echo "MATRIX DONE"
