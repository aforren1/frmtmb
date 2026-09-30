#!/usr/bin/env bash
# Reviewer's driver: every "<pkg> <abs file> <arm> <tag>" line of $1 in
# its own Rscript process, $2 at a time; logs in
# dev/aterms2-rev-suite-<tag>/<arm>__<file>.txt. Env: $3 = gated (0/1),
# $4 = record file prefix for FRMTMB_BRMSPORT_RECORD (optional).
WT=/c/Users/adf44/source/r/frmtmb-wt-aterms2
LIST=$1
P=${2:-12}
GATED=${3:-0}
REC=${4:-}
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export TMP="C:\\Users\\adf44\\AppData\\Local\\Temp\\1"
export TEMP="$TMP"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/aterms2-rev-stan-cache"
export NOT_CRAN=true
if [ "$GATED" = "1" ]; then export FRMTMB_BRMS_FIT_TESTS=true; fi
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
while read -r pkg f arm tag; do
  [ -z "$pkg" ] && continue
  out=$WT/dev/aterms2-rev-suite-$tag
  mkdir -p "$out"
  while [ "$(jobs -rp | wc -l)" -ge "$P" ]; do sleep 1; done
  log="$out/${arm}__$(basename "$f" .R).txt"
  (
    cd "$(dirname "$f")" || exit 1
    if [ -n "$REC" ]; then
      export FRMTMB_BRMSPORT_RECORD="${REC}-${tag}-${arm}-$(basename "$f" .R).tsv"
    fi
    "$R" "$WT/dev/aterms2-rev-testfile.R" "$pkg" "$f" "$arm" > "$log" 2>&1
  ) &
done < "$LIST"
wait
echo "done $(grep -l '^RESULT' $WT/dev/aterms2-rev-suite-*/*.txt | wc -l) logs with RESULT"
