#!/usr/bin/env bash
# Reviewer of lane seflags (check of punch 2b): the core, frmtmb.sample
# and frmtmb.spline test suites on the lane library, one R process per
# file, every gate on (as dev/ciharden-scan.sh sets them), through
# dev/ciharden-run1.R. Logs: dev/seflags-rev4-log/suite-<cfg>/.
#   bash dev/seflags-rev4-suite.sh <ref|ob> <N>
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
cfg="$1"; N="${2:-7}"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$W/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
unset FRMTMB_SAMPLER_GATES
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
       unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS ;;
  ob) RS="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
      export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 ;;
esac
LIB=C:/Users/adf44/source/r/wt-seflags-lib
OUT=dev/seflags-rev4-log/suite-$cfg
rm -rf "$OUT"; mkdir -p "$OUT"
list=$OUT.jobs; : > "$list"
for f in tests/testthat/test-*.R; do echo "frmtmb $f" >> "$list"; done
for e in sample spline; do
  for f in extensions/frmtmb.$e/tests/testthat/test-*.R; do
    echo "frmtmb.$e $f" >> "$list"
  done
done
order=$OUT.order
# longest first, by the lane's ref8 scan (read only)
for j in $(seq 1 $(wc -l < "$list")); do
  read -r p f < <(sed -n "${j}p" "$list")
  t=$(grep -a "^ELAPSED" "dev/ciharden-log/seflags-ref8-suite/$p--$(basename "$f" .R).txt" 2>/dev/null | awk '{print $2}')
  echo "${t:-999999} $p $f"
done | sort -k1,1nr | cut -d" " -f2- > "$order"
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  ( "$RS" dev/ciharden-run1.R "$LIB" "$p" "$f" < /dev/null > "$OUT/$p--$(basename "$f" .R).txt" 2>&1 ) &
done < "$order"
wait
grep -h "^RESULT" "$OUT"/*.txt > "$OUT.sum"
echo "ran $(ls "$OUT"/*.txt | wc -l) of $(wc -l < "$list")" >> "$OUT.sum"
