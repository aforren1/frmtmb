#!/usr/bin/env bash
# Short check, item 8: dev/predread-audit2.R on the review's own mutant builds of the final
# tree, in rev-predread-lib, against the unmutated build. Three at a time.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
O=dev/review-predread-out/r2/mut
mkdir -p "$O"
sed 's#dev/predread-log/audit2-%s.tsv#dev/review-predread-out/r2/mut/audit2-%s.tsv#; s#dev/predread-log/audit2-#dev/review-predread-out/r2/mut/audit2-#' dev/predread-audit2.R > dev/review-predread-audit2-copy.R
"$RS" dev/review-predread-audit2-copy.R C:/Users/adf44/source/r/wt-predread-lib base > "$O/base.txt" 2>&1
one() {
  nm="$1"
  "$RS" dev/review-predread-mutants.R "$nm" > "$O/$nm-build.txt" 2>&1
  "$RS" dev/review-predread-audit2-copy.R "C:/Users/adf44/source/r/rev-predread-lib/$nm" "$nm" > "$O/$nm.txt" 2>&1
}
for nm in cs_se cs_se_rev cs_factor_nd me_var_double me_var_drop me_half_sd me_shared_z me_gr_perrow me_store_off emm_me_design link_identity curve_check_ratio resp_read_off emm_cs_design; do
  one "$nm" &
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 2; done
done
wait
echo finished > "$O/done"
