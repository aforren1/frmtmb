#!/usr/bin/env bash
# Reviewer of lane optima (re-check): every test file of core and the seven
# extensions, one R process per file, N at a time, through
# dev/optima-rev-run-tests.R (the lane's core, nlminb_best_par() traced).
#   bash dev/optima-rev2-run-par.sh <tier> [N] [extensions, default all]
# Output under dev/optima-rev2-out/<tier>/ (gitignored): one .txt per
# file, events.tsv per file, and <tier>.log with the RESULT lines.
set -u
tier="$1"; N="${2:-16}"; EXT="${3:-coupling eam latent learn ode sample spline}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="${OPTREV_RSCRIPT:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
RUN="$ROOT/dev/optima-rev-run-tests.R"
OUT="$ROOT/dev/optima-rev2-out/$tier"
LOG="$ROOT/dev/optima-rev2-out/$tier.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/optima-rev2-out/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
list="$OUT/jobs.lst"
for f in "$ROOT"/tests/testthat/test-*.R; do echo "frmtmb $f"; done > "$list"
for p in $EXT; do
  for f in "$ROOT"/extensions/frmtmb.$p/tests/testthat/test-*.R; do
    echo "frmtmb.$p $f"
  done
done >> "$list"
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$p--$(basename "$f" .R)"
  ( "$RS" "$RUN" "$p" "$f" "$OUT/$b.events" > "$OUT/$b.txt" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  want=$((want + 1))
  b="$p--$(basename "$f" .R)"
  hit=$(grep -a "^RESULT" "$OUT/$b.txt" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else
    echo "$p NO RESULT LINE for $f" >> "$LOG"
  fi
done < "$list"
echo "ran $ran of $want" >> "$LOG"
echo "ran $ran of $want"
