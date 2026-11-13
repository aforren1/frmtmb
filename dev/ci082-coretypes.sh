#!/usr/bin/env bash
# Run dev/ci082-probe.R for one case over OpenBLAS builds and core types,
# at most 3 R processes at once.
#   bash dev/ci082-coretypes.sh <lib|base> <case> <k> <ver...>
# Output: dev/ci082-out/<arm>-<case>-k<k>-<ver>-<coretype>.txt
set -u
ARM="$1"; CASE="$2"; K="$3"; shift 3
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
mkdir -p dev/ci082-out
jobs_run=0
for VER in "$@"; do
  for CT in default Haswell SkylakeX Zen Sandybridge; do
    if [ "$VER" = ref ]; then
      R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
      [ "$CT" = default ] || continue
    else
      R="dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"
    fi
    OUT="dev/ci082-out/$ARM-$CASE-k$K-$VER-$CT.txt"
    (
      if [ "$CT" = default ]; then unset OPENBLAS_CORETYPE; else
        export OPENBLAS_CORETYPE="$CT"; fi
      "$R" dev/ci082-probe.R "$ARM" "$CASE" "$K" > "$OUT" 2>&1
      echo "exit $?" >> "$OUT"
    ) &
    jobs_run=$((jobs_run + 1))
    if [ "$jobs_run" -ge "${MAXJ:-3}" ]; then wait -n; jobs_run=$((jobs_run - 1)); fi
  done
done
wait
grep -H "RESULT\|^exit" dev/ci082-out/$ARM-$CASE-k$K-*.txt
