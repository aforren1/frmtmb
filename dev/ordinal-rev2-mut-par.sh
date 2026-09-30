#!/usr/bin/env bash
# Reviewer re-check (punch round 1), lane ordinal: "<arm> <pkg> <file>" jobs, one R process each,
# N at a time, through dev/ordinal-rev-mutrun.R.
# Output: dev/ordinal-rev2-suite-mut/, mut.log with the RESULT lines.
set -u
list="$1"; N="${2:-20}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordinal"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/ordinal-rev-mutrun.R"
OUT="$ROOT/dev/ordinal-rev2-suite-mut"
LOG="$OUT/mut.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r a p f; do
  [ -n "$a" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  tag=$(basename "$a")
  o="$OUT/$tag--$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" "$a" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r a p f; do
  [ -n "$a" ] || continue
  want=$((want + 1)); tag=$(basename "$a")
  o="$OUT/$tag--$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then ran=$((ran + 1)); echo "$tag $hit" >> "$LOG"
  else echo "$tag NO RESULT LINE for $p $f" >> "$LOG"; fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
echo "RAN $ran of $want"
