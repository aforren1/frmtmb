#!/usr/bin/env bash
# Short check of lane ci082's fix round: test files on the lane library,
# reference BLAS and OpenBLAS 0.3.32 + LAPACK (Sandybridge), at most
# JOBS (default 2) R processes at once. test-se-limit.R runs 1 thread.
#   bash dev/reviews/ci082-short-matrix.sh <outdir> <cfgs> <file...>
set -u
OUT="$1"; CFGS="$2"; shift 2
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export NOT_CRAN=true
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
    case "$F" in *test-se-limit.R) TH=1;; *) TH=4;; esac
    LOG="$OUT/$(basename "$F" .R)-$CFG.txt"
    (
      export OPENBLAS_NUM_THREADS="$TH"
      if [ "$CT" = default ]; then unset OPENBLAS_CORETYPE; else
        export OPENBLAS_CORETYPE="$CT"; fi
      echo "threads: $TH" > "$LOG"
      "$R" dev/ci082-run1.R lib frmtmb "$F" >> "$LOG" 2>&1
      echo "exit $?" >> "$LOG"
    ) &
    n=$((n + 1))
    if [ "$n" -ge "${JOBS:-2}" ]; then wait -n; n=$((n - 1)); fi
  done
done
wait
grep -H "RESULT\|^exit\|^lib:\|^R:\|^threads:" "$OUT"/*.txt
echo "MATRIX DONE"
