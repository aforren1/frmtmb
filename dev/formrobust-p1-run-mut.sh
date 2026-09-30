#!/usr/bin/env bash
# Punch round 1 mutation driver: one R process per (mutant, test file),
# through the reviewer's runner, plus each test file unmutated.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-formrobust
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/formrobust-rev-run-tests.R"
LIST="$ROOT/dev/formrobust-log/p1-jobs-mut.txt"
OUT="$ROOT/dev/formrobust-log/p1-mut-files"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r m p f; do
  [ -n "$m" ] || continue
  mp="$ROOT/dev/formrobust-rev-mut/$m.R"
  [ -f "$mp" ] || mp="$ROOT/dev/formrobust-p1-mut/$m.R"
  ( REVMUTANT="$mp" "$RS" "$RUN" "$p" "$f" > "$OUT/$m--$(basename "$f" .R).txt" 2>&1 ) &
done < "$LIST"
cut -d' ' -f2,3 "$LIST" | sort -u | while read -r p f; do
  ( "$RS" "$RUN" "$p" "$f" > "$OUT/CONTROL--$(basename "$f" .R).txt" 2>&1 ) &
done
wait
for o in "$OUT"/*.txt; do
  echo "$(basename "$o" .txt) $(grep -a '^MUTATED\|^Error' "$o" | head -2 | tr '\n' ' ') $(grep -a '^RESULT' "$o" | tail -1)"
done > "$ROOT/dev/formrobust-log/p1-mut.log"
cat "$ROOT/dev/formrobust-log/p1-mut.log"
