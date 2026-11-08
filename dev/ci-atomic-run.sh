#!/usr/bin/env bash
# One test file per R process, 3 at a time, on the lane library, for the
# test files that can reach stan_cache_write() (dev/ci-cachecallers.py).
#   bash dev/ci-atomic-run.sh <job list> <out dir>
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-ci
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE="$(pwd -W)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
OUT=$2; mkdir -p "$OUT"
L="C:/Users/adf44/source/r/wt-ci-lib;C:/Users/adf44/source/r/rellib-r18"
n=0
while read -r pkg f; do
  n=$((n + 1))
  echo "$n $pkg $f"
done < "$1" | xargs -P 3 -L 1 bash -c '
  R_PARALLEL_PORT=$((11900 + $0)) "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" \
    dev/ci-runfile.R "$2" "$1" "'"$L"'" > "'"$OUT"'/$0-$(basename $2 .R).log" 2>&1'
grep -h RESULT "$OUT"/*.log | sort
