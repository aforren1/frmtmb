#!/usr/bin/env bash
# Records the names of every fit of every test file of one package on
# one arm, one R process per file, N at a time
# (dev/brmsnames-1024-suitenames.R).
#
#   bash dev/brmsnames-1024-suitenames.sh base|lane <pkg> <test dir> [N]
#
# Log: dev/brmsnames-1024-log/suitenames-<arm>-<pkg>.log, one RESULT
# line per file, written by each process's own redirection.
set -u
arm="$1"; pkg="$2"; dir="$3"; N="${4:-8}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
OUT="$ROOT/dev/brmsnames-1024-log/suitenames-$arm-$pkg"
mkdir -p "$OUT"
for f in "$dir"/test-*.R; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$(basename "$f" .R)"
  "$RS" "$ROOT/dev/brmsnames-1024-suitenames.R" "$arm" "$pkg" "$f" \
    > "$OUT/$b.txt" 2>&1 &
done
wait
grep -h "^RESULT" "$OUT"/*.txt > "$ROOT/dev/brmsnames-1024-log/suitenames-$arm-$pkg.log"
echo "files: $(ls "$OUT"/*.txt | wc -l) results: $(wc -l < "$ROOT/dev/brmsnames-1024-log/suitenames-$arm-$pkg.log")"
