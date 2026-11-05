#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/review-latres-out/mut
for m in none ma_flip arma_flip start_var cosy_swap condar_sign no_on_response all_on_response cdf_ignored no_guard methods_no_shift new_rows_no_draw refuse_off sderr_scale2 row_order; do
  for t in test-latres test-latres-brms; do
    echo "$m tests/testthat/$t.R dev/review-latres-out/mut/$m--$t.log"
  done
done | xargs -P 4 -L 1 sh -c '"$0" dev/review-latres-mutants.R "$1" "$2" > "$3" 2>&1' "$R"
grep -h "^RESULT\|^  failed\|^Error" dev/review-latres-out/mut/*.log > dev/review-latres-out/mutants.txt
