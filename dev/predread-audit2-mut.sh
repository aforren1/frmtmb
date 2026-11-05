#!/bin/bash
# Show that dev/predread-audit2.R catches the review's mutants: build each one
# from dev/review-predread-mutants.R's list into the lane's own library tree, run
# the cells, count the cells that differ from the unmutated build.
#   [MUTS="a b"] bash dev/predread-audit2-mut.sh <scratch dir>
# Builds go to wt-predread-lib/mut/<name> and a source copy to the session
# scratchpad, removed after; the reviewer's own directories are not touched.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export R_LIBS="C:/Users/adf44/source/r/wt-predread-lib;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
SRC_ROOT="$1"
LIB_ROOT=C:/Users/adf44/source/r/wt-predread-lib/mut
mkdir -p "$LIB_ROOT" "$SRC_ROOT"
sed "s#^SRC <- .*#SRC <- \"$SRC_ROOT\"#; s#^LIB <- .*#LIB <- \"$LIB_ROOT\"#" \
  dev/review-predread-mutants.R > dev/predread-log/mutants-copy.R
one() {
  nm="$1"
  "$RS" dev/predread-log/mutants-copy.R "$nm" > "dev/predread-log/mut-$nm-build.txt" 2>&1
  "$RS" dev/predread-audit2.R "$LIB_ROOT/$nm" "mut-$nm" \
    > "dev/predread-log/mut-$nm-audit2.txt" 2>&1
}
ALL="cs_se cs_se_rev cs_factor_nd me_var_double me_var_drop me_half_sd me_shared_z me_gr_perrow me_store_off emm_me_design link_identity curve_check_ratio resp_read_off emm_cs_design"
for nm in ${MUTS:-$ALL}; do
  one "$nm" &
  while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 2; done
done
wait
rm -rf "$SRC_ROOT" "$LIB_ROOT"
echo mutants done
