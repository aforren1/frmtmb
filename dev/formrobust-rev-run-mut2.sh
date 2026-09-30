#!/usr/bin/env bash
# Reviewer mutation driver: one R process per (mutant, test file).
# Also runs every test file once unmutated as the control.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-formrobust
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/formrobust-rev-run-tests.R"
LIST="$ROOT/dev/formrobust-rev-log/jobs-mut2.txt"
OUT="$ROOT/dev/formrobust-rev-log/mut2-files"
N="${1:-12}"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/formrobust-rev-stan-cache"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r m p f; do
  [ -n "$m" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$m--$(basename "$f" .R).txt"
  ( REVMUTANT="$ROOT/dev/formrobust-rev-mut/$m.R" "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$LIST"
# controls: each distinct test file unmutated
cut -d' ' -f2,3 "$LIST" | sort -u | while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/CONTROL--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done
wait
for o in "$OUT"/*.txt; do
  echo "$(basename "$o" .txt) $(grep -a '^MUTATED\|^Error' "$o" | head -2 | tr '\n' ' ') $(grep -a '^RESULT' "$o" | tail -1)"
done > "$ROOT/dev/formrobust-rev-log/mut2.log"
echo done
