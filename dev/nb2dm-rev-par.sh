#!/usr/bin/env bash
# Run "<package> <test file>" jobs one R process per file, N at a time,
# through dev/nb2dm-run-test.R (reviewer copy of dev/nb2dm-par.sh), and gather the RESULT lines.
#   bash dev/nb2dm-par.sh <tag> <job list> [N] [lib]
# Per-file logs: dev/nb2dm-rev-out/<tag>-files/; summary: dev/nb2dm-rev-out/<tag>.txt.
set -u
tag="$1"; list="$2"; N="${3:-6}"; lib="${4:-}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="${RSCRIPT:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
OUT="$ROOT/dev/nb2dm-rev-out/$tag-files"
LOG="$ROOT/dev/nb2dm-rev-out/$tag.txt"
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
cd "$ROOT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  if [ -n "$lib" ]; then
    ( "$RS" dev/nb2dm-run-test.R "$p" "$f" "$lib" > "$o" 2>&1 ) &
  else
    ( "$RS" dev/nb2dm-run-test.R "$p" "$f" > "$o" 2>&1 ) &
  fi
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$hit" >> "$LOG"
    grep -a -A20 "^FAILING BLOCK" "$o" >> "$LOG"
  else
    echo "NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
