#!/usr/bin/env bash
# Review of lane boundfam: run a job list ("<pkg> <test file>" lines)
# through an R script, N R processes at a time, one file per process.
#   bash dev/reviews/boundfam-review-par.sh <script.R> <job list> <outdir> [N] [port]
set -u
script="$1"; list="$2"; OUT="$3"; N="${4:-2}"; port="${5:-23100}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/stan-cache"
export NOT_CRAN=true
mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; continue; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  b="$OUT/$p--$(basename "$f" .R)"
  port=$((port + 1))
  ( R_PARALLEL_PORT=$port "$RS" "$script" "$p" "$f" "$b.prof" > "$b.txt" 2>&1 ) &
done < "$list"
wait
echo "done"
