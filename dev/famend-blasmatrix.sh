#!/usr/bin/env bash
# Lane famend, fix round: the perturbation checks and every changed test
# file under reference BLAS and OpenBLAS 0.3.26 / 0.3.32 LAPACK (core
# type unset and Haswell), three R processes at a time.
#   bash dev/famend-blasmatrix.sh
WT=C:/Users/adf44/source/r/frmtmb-wt-famend
cd "$WT"
LIB=C:/Users/adf44/source/r/wt-famend-lib
OUT=dev/famend-out/blas
mkdir -p $OUT
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
files="test-family-fitend.R test-families-extreme.R test-osa-cdf.R
test-cens-trunc.R test-linpred-layers.R test-prior-table-brms.R"
jobs_list=()
for arm in ref ob026 ob026h ob032 ob032h; do
  case $arm in
    ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"; CT="";;
    ob026) RS=dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe; CT="";;
    ob026h) RS=dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe; CT=Haswell;;
    ob032) RS=dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe; CT="";;
    ob032h) RS=dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe; CT=Haswell;;
  esac
  jobs_list+=("$arm|$RS|$CT|perturb")
  for f in $files; do jobs_list+=("$arm|$RS|$CT|$f"); done
done
for j in "${jobs_list[@]}"; do
  IFS='|' read -r arm RS CT what <<< "$j"
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 2; done
  (
    if [ -n "$CT" ]; then export OPENBLAS_CORETYPE=$CT; else unset OPENBLAS_CORETYPE; fi
    if [ "$what" = perturb ]; then
      "$RS" dev/famend-perturb.R $LIB 1:20 > $OUT/$arm-perturb.txt 2>&1
    else
      "$RS" dev/famend-runtest.R $LIB frmtmb tests/testthat/$what \
        > $OUT/$arm-$what.txt 2>&1
    fi
  ) &
done
wait
for f in $OUT/*.txt; do
  echo "== $(basename $f)"; grep -a "^RESULT\|^PERTURB\|FAIL\|expectation_" $f
done
