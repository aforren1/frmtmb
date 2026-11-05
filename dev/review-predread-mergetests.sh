#!/usr/bin/env bash
# Item 8: the tests the three lanes changed or added, on the trial merge, one per process.
set -u
cd /c/Users/adf44/source/r/rev-predread-src/merge
P=/c/Users/adf44/source/r/frmtmb-wt-predread
O=$P/dev/review-predread-out/merge-tests
mkdir -p "$O"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_STAN_CACHE=$P/dev/stan-cache
LIB=C:/Users/adf44/source/r/rev-predread-lib/merge
N=${1:-6}
while read -r f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  case "$f" in
    extensions/frmtmb.sample/*) p=frmtmb.sample ;;
    *) p=frmtmb ;;
  esac
  ( "$RS" "$P/dev/review-predread-runone.R" "$p" "$f" "$LIB" > "$O/$p--$(basename "$f" .R).txt" 2>&1 ) &
done < ../merge.tests
wait
grep -ah "^RESULT\|^FAILED" "$O"/*.txt | cut -c1-260 > "$O/../merge-tests.sum"
echo finished >> "$O/../merge-tests.sum"
