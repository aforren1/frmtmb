#!/usr/bin/env bash
# The record of dev/rel074-saved-mixed.R, rerun on this lane's core: a
# saved cox() and car() fit read in four libraries that hold this
# lane's frmtmb and neither extension, frmtmb.survival only,
# frmtmb.spatial only, or both (the extensions are rellib-r12's 0.74.0
# copies; this lane changes no R code of theirs).
#   bash dev/smallfix-saved-mixed.sh <arm: before|after>
# before reads with rellib-r12's frmtmb, after with the lane library's.
# Log: dev/smallfix-out/saved-mixed-<arm>.txt
set -u
arm=$1
WT=C:/Users/adf44/source/r/frmtmb-wt-smallfix
R12=C:/Users/adf44/source/r/rellib-r12
LANE=C:/Users/adf44/source/r/wt-smallfix-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
core=$LANE; [ "$arm" = before ] && core=$R12
L=$WT/dev/smallfix-out/saved-mixed-libs-$arm
rm -rf "$L"; mkdir -p "$L"
for s in corelib lib-surv lib-spat lib-both; do
  mkdir -p "$L/$s"; cp -r "$core/frmtmb" "$L/$s/"
done
for s in lib-surv lib-both; do cp -r "$R12/frmtmb.survival" "$L/$s/"; done
for s in lib-spat lib-both; do cp -r "$R12/frmtmb.spatial" "$L/$s/"; done
out=$WT/dev/smallfix-out/saved-mixed-$arm.txt
rds=$WT/dev/smallfix-out/saved-mixed.rds
# the fit is saved once, with rellib-r12, as the 0.74.0 record saved it
[ -f "$rds" ] || "$RS" "$WT/dev/rel074-saved-mixed.R" save "$rds" > "$out.save" 2>&1
cat "$out.save" > "$out"
for s in corelib lib-surv lib-spat lib-both; do
  echo "== read with $s" >> "$out"
  "$RS" "$WT/dev/rel074-saved-mixed.R" read "$rds" "$L/$s" >> "$out" 2>&1
done
rm -rf "$L"
cat "$out"
