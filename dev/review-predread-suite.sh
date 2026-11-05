#!/usr/bin/env bash
# One file per process over a job list, N at a time, with the lane's library.
#   bash dev/review-predread-suite.sh <tag> <joblist> <N>
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
tag="$1"; list="$2"; N="$3"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="dev/review-predread-out/$tag-files"
mkdir -p "$OUT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-predread/dev/stan-cache
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB="${FRMTMB_RELLIB:-C:/Users/adf44/source/r/wt-predread-lib}"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  [ -s "$o" ] && grep -qa "^RESULT" "$o" && continue
  ( "$RS" dev/release/run-tests.R "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
echo done
