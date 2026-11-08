#!/usr/bin/env bash
# Review of lane arstep: second chain, one R process at a time.
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
W=C:/Users/adf44/source/r/frmtmb-wt-arstep
O=$W/dev/review-arstep-out
cd "$W" || exit 1
"$R" dev/review-arstep-mut.R nomask > "$O/mut-nomask.txt" 2>&1
for m in control nosettle norestore; do
  "$R" dev/review-arstep-mut.R "$m" "$W/dev/review-arstep-settle-test.R" > "$O/settle-$m.txt" 2>&1
done
"$R" dev/review-arstep-restore.R > "$O/restore.txt" 2>&1
"$R" dev/review-arstep-blas.R refblas > "$O/blas-ref.txt" 2>&1
OPENBLAS_NUM_THREADS=4 "$W/dev/review-arstep-ob/Rob0.3.26-lapack/bin/x64/Rscript.exe" dev/review-arstep-blas.R ob0.3.26 > "$O/blas-ob26.txt" 2>&1
OPENBLAS_NUM_THREADS=4 "$W/dev/review-arstep-ob/Rob0.3.32-lapack/bin/x64/Rscript.exe" dev/review-arstep-blas.R ob0.3.32 > "$O/blas-ob32.txt" 2>&1
R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win "$R" dev/review-arstep-brmsfit.R > "$O/brmsfit.txt" 2>&1
"$R" dev/review-arstep-cost.R a300,p600,p1000,pois120,pois240,bern240 2 > "$O/cost-1.txt" 2>&1
echo CHAIN2DONE > "$O/chain2-done.txt"
