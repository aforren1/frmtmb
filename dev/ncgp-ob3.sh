#!/usr/bin/env bash
# Row 10a's gradient on the reference BLAS and on OpenBLAS with LAPACK, 3 runs each.
O=dev/ncgp-out/ob-row10a.txt
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/ncgp-out/ob26-build.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/ncgp-out/ob32-build.txt 2>&1
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
: > "$O"
for rep in 1 2 3; do
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/ncgp-row10a.R >> "$O" 2>&1
  for th in 1 4; do
    for v in 26 32; do
      echo "== ob$v threads $th rep $rep" >> "$O"
      OPENBLAS_NUM_THREADS=$th dev/ciharden-out/Rob0.3.$v-lapack/bin/Rscript.exe dev/ncgp-row10a.R >> "$O" 2>&1
    done
  done
done
echo OB3-DONE >> "$O"
