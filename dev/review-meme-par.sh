#!/usr/bin/env bash
# Review of lane meme: run "<pkg> <file> <arm>" jobs N at a time, one R
# process per file. Usage: bash dev/review-meme-par.sh <jobs> <tag> [N]
set -u
list="$1"; tag="$2"; N="${3:-6}"
WT=/c/Users/adf44/source/r/frmtmb-wt-meme
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$WT/dev/review-meme-out/$tag"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
mkdir -p "$OUT"
while read -r p f arm; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  o="$OUT/$arm--$p--$(basename "$f" .R).txt"
  ( cd "$WT" && "$RS" "$WT/dev/review-meme-run1.R" "$p" "$f" "$arm" > "$o" 2>&1 ) &
done < "$list"
wait
grep -a -h "^RESULT" "$OUT"/*.txt
for o in "$OUT"/*.txt; do grep -q "^RESULT" "$o" || echo "NO RESULT: $o"; done
echo "DONE $tag"
