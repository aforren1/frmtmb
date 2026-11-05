#!/usr/bin/env bash
# bash dev/nlfix-suite.sh <tag> <gated 0|1> <N> <pkg> [<pkg> ...]
# Every test file of each package, one R process per file, N at once, on the
# lane's library. Writes dev/nlfix-suite-<tag>/<pkg>--<file>.txt and prints a
# table of files that are not clean. Rtools at the END of PATH.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
tag=$1; gated=$2; N=$3; shift 3
out=dev/nlfix-suite-$tag
mkdir -p "$out"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-nlfix/dev/stan-cache"
if [ "$gated" = 1 ]; then export FRMTMB_BRMS_FIT_TESTS=true; fi
: > "$out/jobs.txt"
for p in "$@"; do
  if [ "$p" = frmtmb ]; then d=tests/testthat; else d=extensions/$p/tests/testthat; fi
  for f in "$d"/test-*.R; do
    if [ -n "${ONLY:-}" ] && ! echo "$f" | grep -qE "$ONLY"; then continue; fi
    echo "$p $f" >> "$out/jobs.txt"
  done
done
cat "$out/jobs.txt" | xargs -P "$N" -L 1 bash -c \
  'p=$0; f=$1; export R_PARALLEL_PORT=$((12200 + RANDOM % 1700)); b=$(basename "$f" .R); o="'"$out"'/$p--$b"; if [ ! -s "$o.txt" ]; then timeout 3000 "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/nlfix-testrun.R "$p" "$f" "$o.txt" > "$o.err" 2>&1 || echo "PROCESS FAILED $p $f" >> "$o.txt"; fi'
grep -h "^FILE\|PROCESS FAILED" "$out"/*.txt | awk '{ n++; split($3,a,"="); split($4,b,"="); split($5,c,"="); split($6,d,"="); split($7,e,"="); pa+=a[2]; fa+=b[2]; er+=c[2]; sk+=d[2]; wa+=e[2]; if (b[2]>0 || c[2]>0 || e[2]>0 || $1=="PROCESS") print } END { printf "TOTAL files=%d PASS=%d FAIL=%d ERROR=%d SKIP=%d WARN=%d\n", n, pa, fa, er, sk, wa }'
