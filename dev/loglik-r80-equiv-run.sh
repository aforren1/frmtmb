#!/usr/bin/env bash
# Remeasure the sampled pairs that carry a prior (dev/loglik-r80-equiv-ids.tsv)
# on the lane build, three at a time, into dev/loglik-equiv-out.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-loglik
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-loglik-lib;C:/Users/adf44/source/r/rellib-r18"
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/stan-cache"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=true
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=dev/loglik-equiv-out
mkdir -p "$OUT"
N=3
# the jobs must be children of this shell, or `wait` returns at once
while read -r id; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  ( "$RS" dev/equiv-measure.R "$id" "$OUT" > "$OUT/$id.log" 2>&1 ) &
done < <(tail -n +2 dev/loglik-r80-equiv-ids.tsv | cut -f1)
wait
echo "equiv run done"
