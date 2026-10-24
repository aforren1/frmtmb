#!/usr/bin/env bash
# Each lane's fit and curve comparisons, rerun on the merged build: the
# "before" arm on rellib-r11 (frmtmb 0.73.0 and its seven extensions),
# the "after" arm on rellib-r12 (frmtmb 0.74.0 and its eight). The
# lanes' own scripts are run, unedited where they take their libraries
# as arguments (mvnlatent, spatial). splitsurv's dumps name the lane
# library and write to dev/splitsurv-out/, so they run from copies in
# dev/rel074-log/eq/ with those two strings replaced, which sed prints.
# Every arm is one process, so a difference is a difference of builds
# (dev/rtmb-pitfalls.md item 21).
#
#   bash dev/rel074-equality.sh > dev/rel074-log/equality.txt 2>&1
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R11=C:/Users/adf44/source/r/rellib-r11
R12=C:/Users/adf44/source/r/rellib-r12
E=dev/rel074-log/eq
rm -rf "$E"; mkdir -p "$E/out"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache

echo "== splitsurv: curves, fits (dev/splitsurv-{curve,fit}-dump.R)"
for f in splitsurv-curve-dump.R splitsurv-fit-dump.R splitsurv-compare.R \
         splitsurv-fit-compare.R; do
  sed -e "s#C:/Users/adf44/source/r/wt-splitsurv-lib#$R12#g" \
      -e "s#dev/splitsurv-out#$E/out#g" "dev/$f" > "$E/$f"
  diff "dev/$f" "$E/$f" | grep '^>' | sed "s/^/   $f: /"
done
# the curve compare also reads the compat dumps; the merged table has
# its own check (dev/rel074-compat.R), so a copy of the lane's two
# dumps stands in
cp ../frmtmb-wt-splitsurv/dev/splitsurv-out/compat-before.rds \
   ../frmtmb-wt-splitsurv/dev/splitsurv-out/compat-after.rds "$E/out/"
( "$RS" "$E/splitsurv-curve-dump.R" before "$E/out/curve-before.rds" &&
  "$RS" "$E/splitsurv-fit-dump.R" before "$E/out/fit-before.rds" ) \
  > "$E/ss-before.log" 2>&1 &
( "$RS" "$E/splitsurv-curve-dump.R" after "$E/out/curve-after.rds" &&
  "$RS" "$E/splitsurv-fit-dump.R" after "$E/out/fit-after.rds" ) \
  > "$E/ss-after.log" 2>&1 &
echo "== mvnlatent: 59 mixture_mvn() fits (dev/mvnlatent-fitcompare.R)"
"$RS" dev/mvnlatent-fitcompare.R base "$E/out/mvn-base.rds" \
  > "$E/mvn-base.log" 2>&1 &
"$RS" dev/mvnlatent-fitcompare.R "$R12" "$E/out/mvn-merged.rds" \
  > "$E/mvn-merged.log" 2>&1 &
echo "== spatial: car() and spde() fits, frm_sample() draws"
( SPATIAL_LIBS="$R11" SPATIAL_PKGS="" \
    "$RS" dev/spatial-fitcompare.R fit "$E/out/sp-before.rds" &&
  SPATIAL_LIBS="$R11" SPATIAL_PKGS="" \
    "$RS" dev/spatial-sample-check.R "$E/out/sps-before.rds" ) \
  > "$E/sp-before.log" 2>&1 &
( SPATIAL_LIBS="$R12" SPATIAL_PKGS="frmtmb.spatial" \
    "$RS" dev/spatial-fitcompare.R fit "$E/out/sp-after.rds" &&
  SPATIAL_LIBS="$R12" SPATIAL_PKGS="frmtmb.spatial" \
    "$RS" dev/spatial-sample-check.R "$E/out/sps-after.rds" ) \
  > "$E/sp-after.log" 2>&1 &
wait
for l in "$E"/*.log; do echo "-- $l"; grep -E "^lib|lib:|frmtmb [0-9]" "$l" | head -6; done
echo "== splitsurv compare"
"$RS" "$E/splitsurv-compare.R" 2>&1 | sed -n '/== curve outputs/,$p' | grep -v "^== frm_compat"
"$RS" "$E/splitsurv-fit-compare.R" 2>&1
echo "== mvnlatent compare"
"$RS" dev/mvnlatent-fitcompare-diff.R "$E/out/mvn-base.rds" "$E/out/mvn-merged.rds" 2>&1 | head -20
echo "== spatial compare"
"$RS" dev/spatial-fitcompare.R compare "$E/out/sp-before.rds" "$E/out/sp-after.rds" 2>&1
"$RS" -e 'a <- readRDS("'"$E"'/out/sps-before.rds"); b <- readRDS("'"$E"'/out/sps-after.rds"); for (k in intersect(names(a), names(b))) { x <- a[[k]]; y <- b[[k]]; cat(k, "identical:", identical(x, y), "\n") }' 2>&1
echo "EQUALITY DONE"
