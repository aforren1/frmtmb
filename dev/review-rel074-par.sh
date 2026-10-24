#!/usr/bin/env bash
# Release review of 0.74.0: run a job list ("<package> <test file>" per
# line) one file per R process, N at a time, through
# dev/release/run-tests.R (which reads rellib-r12), writing logs only
# under dev/review-rel074-log/<tier>/ (gitignored).
#   bash dev/review-rel074-par.sh <tier> <job list> [N]
set -u
tier="$1"; list="$2"; N="${3:-4}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/release/run-tests.R"
OUT="$ROOT/dev/review-rel074-log/$tier"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
export NOT_CRAN=true
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
for o in "$OUT"/*.txt; do
  echo "$(basename "$o" .txt): $(grep -a '^lib:' "$o" | head -1 | sed 's/.*frmtmb [0-9.]* //') | $(grep -a '^RESULT' "$o" | tail -1)"
done > "$OUT.sum"
echo "DONE $tier"
