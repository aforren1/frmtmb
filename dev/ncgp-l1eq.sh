#!/usr/bin/env bash
# dev/ncgp-l1eq.R on the survey's five pairs with an exact gp(), spec names as
# dev/equiv-specs.R keys them, one at a time.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
O=dev/ncgp-out/l1eq.txt
: > "$O"
for k in "gp_pred() x mvbf" "gp_pred() x prior" "lognormal x gp" "student x gp" "gp x nl"; do
  "$RS" dev/ncgp-l1eq.R "$k" 2>&1 | grep -a -E 'NCGPE|^Error' >> "$O"
done
echo L1EQ-DONE
