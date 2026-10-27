#!/usr/bin/env bash
# Review, short check of release fixes R4 to R9: the job list
# dev/review-rel076-out/short-jobs.txt ("<package> <file>"), one file per
# process, 4 at a time, every gate on, on rellib-r14; Rtools last on PATH.
#   bash dev/review-rel076-short.sh
REL=/c/Users/adf44/source/r/frmtmb-wt-rel076
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$REL/dev/review-rel076-out/short
mkdir -p "$OUT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$REL/dev/review-rel076-out/stan-cache
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r14
cd "$REL" || exit 1
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 2; done
  ( "$RS" "$REL/dev/ciharden-run1.R" base "$p" "$REL/$f" \
      > "$OUT/$p--$(basename "$f" .R).log" 2>&1 ) &
done < "$REL/dev/review-rel076-out/short-jobs.txt"
wait
for f in "$OUT"/*--*.log; do
  echo "$(basename "$f") | $(grep -a '^RESULT' "$f" | tail -1) | $(grep -a '^lib:' "$f" | head -1)"
done > "$OUT/summary.txt"
echo ALLDONE >> "$OUT/summary.txt"
