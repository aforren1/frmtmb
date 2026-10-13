#!/usr/bin/env bash
# Lane seflags: the item scripts (fixture 1, the flat-peak sd, the alias
# construction, the simplex face, Malingering_2, the degree-7
# polynomial, the ridge spectrum) on the lane build and the base, under
# one R build. Logs: dev/seflags-log/ext-<tag>/.
#   bash dev/seflags-extras.sh <tag> <ref|ob26>
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
tag="$1"; cfg="$2"
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe" ;;
  ob26) RS="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
        export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 ;;
esac
LL=C:/Users/adf44/source/r/wt-seflags-lib
R7=C:/Users/adf44/source/r/rellib-r7
L="$W/dev/seflags-log/ext-$tag"
S="$W/dev/seflags-out/scripts"
mkdir -p "$L"
bash dev/seflags-retarget.sh dev/relrev069-flatsd.R dev/relrev069-dh7diag.R \
  dev/optima-segap.R
# the base arm of the segap script: rellib-r7 alone
sed 's#"C:/Users/adf44/source/r/wt-seflags-lib",##' "$S/optima-segap.R" \
  > "$S/optima-segap-base.R"
for who in lane base; do
  if [ $who = lane ]; then lib=$LL; else lib=$R7; fi
  "$RS" "$S/relrev069-flatsd.R" "$lib" > "$L/flatsd-$who.txt" 2>&1 &
  "$RS" "$S/relrev069-dh7diag.R" "$lib" > "$L/dh7diag-$who.txt" 2>&1 &
  "$RS" dev/seflags-alias.R "$lib" re > "$L/alias-re-$who.txt" 2>&1 &
  "$RS" dev/seflags-alias.R "$lib" > "$L/alias-$who.txt" 2>&1 &
  "$RS" dev/seflags-bcmmal.R "$lib" > "$L/bcmmal-$who.txt" 2>&1 &
  "$RS" dev/seflags-deg7diag.R "$lib" > "$L/deg7diag-$who.txt" 2>&1 &
  "$RS" dev/seflags-spectrum.R "$lib" 1:8 > "$L/spectrum-$who.txt" 2>&1 &
  wait
done
"$RS" "$S/optima-segap.R" > "$L/segap-lane.txt" 2>&1 &
"$RS" "$S/optima-segap-base.R" > "$L/segap-base.txt" 2>&1 &
"$RS" dev/seflags-ridgecurv.R "$LL" > "$L/ridgecurv-lane.txt" 2>&1 &
wait
echo EXTRAS DONE
