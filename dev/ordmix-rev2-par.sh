#!/usr/bin/env bash
# Reviewer of lane ordmix: one R process per test file, N at a time.
#   bash dev/ordmix-rev-par.sh <arm> <tag> <job list> [N]
# Job list lines: "<package> <test file path>". Output:
# dev/ordmix-rev2-suite-<tag>/ with one log per file and <tag>.log.
set -u
arm="$1"; tag="$2"; list="$3"; N="${4:-20}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/ordmix-rev2-runtest.R"
OUT="$ROOT/dev/ordmix-rev2-suite-$tag"
LOG="$OUT/$tag.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
: "${TMP:?TMP is unset}"
export TMPDIR="$TMP"
export R_PROFILE_USER="$ROOT/dev/ordmix-rev2-profile.R"
export ORDMIX_REV2_LOG="$ROOT/dev/ordmix-rev2-cc-log"
mkdir -p "$ORDMIX_REV2_LOG"
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" "$arm" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else
    echo "$p NO RESULT LINE for $f" >> "$LOG"
  fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
echo "RAN $ran of $want"
