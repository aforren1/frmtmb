#!/usr/bin/env bash
# Short check of lane ci082's fix round: the review's vc cases on the lane
# library, and the attack cases on lane and 0.82.0, three R processes at
# most (one sequential chain each), reference BLAS.
#   bash dev/reviews/ci082-short-run.sh
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/reviews/ci082-short-out
mkdir -p "$O"
(
  for c in ouind csstart cslow fix999; do
    "$R" dev/reviews/ci082-review-vc.R lib "$c" > "$O/lib-$c.txt" 2>&1
  done
) &
(
  for c in usclamp oufar fixclamp invpath; do
    "$R" dev/reviews/ci082-short-attack.R lib "$c" > "$O/atk-lib-$c.txt" 2>&1
  done
) &
(
  for c in usclamp oufar fixclamp invpath; do
    "$R" dev/reviews/ci082-short-attack.R base "$c" > "$O/atk-base-$c.txt" 2>&1
  done
) &
wait
echo RUN DONE
