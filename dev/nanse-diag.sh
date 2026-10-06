#!/usr/bin/env bash
# Rerun the test files where the standard-error warning fired, with
# dev/nanse-run-diag.R, N at a time.
#   bash dev/nanse-diag.sh <job list> <out dir> [N]
set -u
LIST="$1"; OUT="$2"; N="${3:-20}"
W=/c/Users/adf44/source/r/frmtmb-wt-nanse
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$W/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true
: "${TMP:?TMP is unset}"
mkdir -p "$W/$OUT"
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$p--$(basename "$f" .R)"
  ( cd "$(dirname "$f")" && "$RS" "$W/dev/nanse-run-diag.R" "$p" "$f" \
      "$W/$OUT/$b.fire" > "$W/$OUT/$b.txt" 2>&1 ) &
done < "$LIST"
wait
echo done > "$W/$OUT/DONE"
