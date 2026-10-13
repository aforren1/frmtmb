#!/usr/bin/env bash
# Reviewer of lane optrap, final check: dev/ciharden-scan.sh with its output under
# dev/optrap-rev3-log/scan-<run>/ so that the lane's own logs are not
# touched. Every gate on, one R process per test file.
#   bash dev/optrap-rev-scan.sh <run> <ref|ob26> <lib dir|base> <tree> [N] [pkgs]
# <tree>: the checkout whose test files run (the lane's worktree, or the
# merged tree dev/optrap-rev-merge). [pkgs]: "core" (default) or "all".
set -u
run="$1"; cfg="$2"; lib="$3"; TREE="$4"; N="${5:-12}"; pkgs="${6:-core}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-optrap
cd "$TREE"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export TMP=/c/Users/adf44/AppData/Local/Temp TEMP=/c/Users/adf44/AppData/Local/Temp
unset FRMTMB_SAMPLER_GATES
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
       unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS ;;
  ob26) RS="$ROOT/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
        export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 ;;
esac
OUT="$ROOT/dev/optrap-rev3-log/scan-$run"
rm -rf "$OUT"; mkdir -p "$OUT"
list="$OUT/jobs"
: > "$list"
for f in tests/testthat/test-*.R; do echo "frmtmb $f" >> "$list"; done
if [ "$pkgs" = all ]; then
  for e in coupling eam latent learn ode sample spline; do
    for f in $ROOT/extensions/frmtmb.$e/tests/testthat/test-*.R; do
      echo "frmtmb.$e $f" >> "$list"
    done
  done
fi
echo "config $cfg Rscript $RS tree $TREE lib $lib" > "$OUT/CONFIG"
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$ROOT/dev/ciharden-run1.R" "$lib" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
SUM="$OUT/SUM"; : > "$SUM"
while read -r p f; do
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -n 1)
  if [ -n "$hit" ] && ! echo "$hit" | grep -q LOADERROR; then
    ran=$((ran + 1)); echo "$p $hit" >> "$SUM"
  else
    echo "$p NO RESULT LINE for $f" >> "$SUM"
  fi
done < "$list"
echo "ran $ran of $want" >> "$SUM"
tail -n 1 "$SUM"
