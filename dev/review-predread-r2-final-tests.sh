#!/usr/bin/env bash
# Short check of fix round 2: on the merge of the finished lanes with F1 to F4, the files
# reaching any of them, the two merged brms suites, and the file nlfix changed last.
set -u
cd /c/Users/adf44/source/r/rev-predread-src/merge
P=/c/Users/adf44/source/r/frmtmb-wt-predread
O=$P/dev/review-predread-out/r2/merge-final
mkdir -p "$O"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/rev-predread-src/stan-cache
LIB=C:/Users/adf44/source/r/rev-predread-lib/merge-final
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 2; done
  ( "$RS" "$P/dev/review-predread-runone.R" "$p" "$f" "$LIB" > "$O/$p--$(basename "$f" .R).txt" 2>&1 ) &
done < "$P/dev/review-predread-out/r2/jobs-final.txt"
wait
echo finished > "$O/../merge-final.done"
