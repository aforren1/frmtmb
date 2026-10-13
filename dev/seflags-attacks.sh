#!/usr/bin/env bash
# Lane seflags: the SE attack set of lanes nanse and setier (and their
# reviews), the singular, separation and short-fit studies and the
# window designs, on the lane build and on the 0.69.0 base, under one R
# build. The scripts are retargeted copies (dev/seflags-retarget.sh).
#
#   bash dev/seflags-attacks.sh <tag> <ref|ob26> [N]
#
# Logs: dev/seflags-log/att-<tag>/<script>-<lane|base>.txt (and .tsv
# for the studies). Never tail a log while it runs.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
tag="$1"; cfg="$2"; N="${3:-12}"
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe" ;;
  ob26) RS="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
        export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 ;;
  *) echo "config: ref or ob26"; exit 1 ;;
esac
: "${TMP:?TMP is unset}"
LL=C:/Users/adf44/source/r/wt-seflags-lib
R7=C:/Users/adf44/source/r/rellib-r7
L="$W/dev/seflags-log/att-$tag"
S="$W/dev/seflags-out/scripts"
mkdir -p "$L"
export SEFLAGS_OUT="$L"
arm_scripts="setier-rev-a-rev-cases:lane:base setier-rev-a-rev-spread2:lane:base
setier-rev-a-rev-grby:lane:base setier-rev-a-rev-pred:lane:base
setier-rev-a-rev2-b1:lane:merge setier-rev-a-rev2-b1-rho:lane:merge
setier-rev-a-rev3-falseloss-sweep:lane:merge
setier-rev-a-rev2-ridge-re:lane:base setier-rev2-window:lane:base
setier-rev-probe:lane:base setier-rev-c0k:lane:base
setier-rev-c0k-sweep:lane:base setier-rev-aexpb:lane:base
setier-rev-deg7:lane:base setier-rev2-bnd:lane:base
setier-rev-smooth2:lane:base setier-rev-fixes-table:merge:release"
lib_scripts="setier-rev2-trap setier-rev2-short setier-rev2-scale2
setier-rev2-item8 setier-rev-rs20 setier-rev3-trapq setier-rev3-quad"
all=""
for x in $arm_scripts; do all="$all dev/${x%%:*}.R"; done
for x in $lib_scripts setier-rev2-sepmiss setier-singular setier-rev-smallsd \
         setier-sep; do
  all="$all dev/$x.R"
done
bash "$W/dev/seflags-retarget.sh" $all
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
# the longest first
for who in lane base; do
  if [ $who = lane ]; then lib=$LL; blib=$LL; else lib=$R7; blib=base; fi
  jobs_wait
  "$RS" "$S/setier-singular.R" "$blib" "$L/singular-$who.tsv" \
    > "$L/singular-$who.txt" 2>&1 &
  jobs_wait
  "$RS" "$S/setier-rev-smallsd.R" "$lib" "$L/smallsd-$who.tsv" 40 \
    > "$L/smallsd-$who.txt" 2>&1 &
  jobs_wait
  "$RS" "$S/setier-sep.R" "$blib" "$L/sep-$who.tsv" > "$L/sep-$who.txt" 2>&1 &
done
for x in $arm_scripts; do
  s=${x%%:*}; r=${x#*:}; la=${r%%:*}; ba=${r#*:}
  extra=""
  [ "$s" = setier-rev-deg7 ] && extra=7
  jobs_wait
  "$RS" "$S/$s.R" "$la" $extra > "$L/$s-lane.txt" 2>&1 &
  jobs_wait
  "$RS" "$S/$s.R" "$ba" $extra > "$L/$s-base.txt" 2>&1 &
done
for s in $lib_scripts; do
  for who in lane base; do
    if [ $who = lane ]; then lib=$LL; else lib=$R7; fi
    jobs_wait
    "$RS" "$S/$s.R" "$lib" > "$L/$s-$who.txt" 2>&1 &
  done
done
for who in lane base; do
  if [ $who = lane ]; then lib=$LL; else lib=$R7; fi
  jobs_wait
  "$RS" "$S/setier-rev2-sepmiss.R" "$lib" budget few \
    > "$L/sepmiss-$who.txt" 2>&1 &
  jobs_wait
  "$RS" "$S/setier-rev2-sepmiss.R" "$lib" bigp > "$L/sepmiss-bigp-$who.txt" 2>&1 &
done
wait
echo ATTACKS DONE
