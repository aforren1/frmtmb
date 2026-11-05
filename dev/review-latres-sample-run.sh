#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/review-latres-out/sample
for j in "pois 3" "bern 2" "binom3 3" "bern_long 3" "pois_ar9 3" "bern_F 2" "pois_F 3"; do
  echo "$j"
done | xargs -P 4 -L 1 sh -c '"$0" dev/review-latres-sample.R "$1" "$2" > "dev/review-latres-out/sample/$1-$2.log" 2>&1' "$R"
