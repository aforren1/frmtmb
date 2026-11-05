#!/usr/bin/env bash
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for m in control no_order no_nugget first_last start_medoid one_dim no_start_fix no_random; do
  while [ "$(jobs -r | wc -l)" -ge 3 ]; do sleep 2; done
  ( "$RS" dev/gpgap-mut2.R "$m" > "dev/gpgap-out/mut2/$m.txt" 2>&1 ) &
done
( "$RS" dev/review-gpgap-lapcost.R "C:/Users/adf44/source/r/wt-gpgap-lib;C:/Users/adf44/source/r/rellib-r16" > dev/gpgap-out/r2/lapcost-now.txt 2>&1 ) &
( "$RS" dev/review-gpgap-lapcost.R "C:/Users/adf44/source/r/rellib-r16" > dev/gpgap-out/r2/lapcost-base.txt 2>&1 ) &
wait
