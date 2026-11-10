#!/usr/bin/env bash
# Reviewer: remeasure a few pairs with the lane's dev/equiv-measure.R on
# the lane build, three R processes at once.
#   bash dev/review-compat-run.sh <outdir> <id> [<id> ...]
set -u
out=$1; shift
cd "$(dirname "$0")/.."
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-compat-lib;C:/Users/adf44/source/r/rellib-r19"
export FRMTMB_STAN_CACHE="$PWD/dev/stan-cache"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$out"
n=0
for id in "$@"; do
  ( "$RS" dev/equiv-measure.R "$id" "$out" --force > "$out/$id.log" 2>&1
    echo "done $id exit $?" ) &
  n=$((n + 1))
  if [ $((n % 3)) -eq 0 ]; then wait; fi
done
wait
grep -h "^pair\|^lib:\|^verdict" "$out"/*.log
