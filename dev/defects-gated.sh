#!/bin/sh
# Lane wt-defects: the files the ungated suite skipped in, rerun with
# the gate open (FRMTMB_BRMS_FIT_TESTS=true), one R process per file, all
# at once. Lists: dev/defects-log/gated-<pkg>-files.txt.
#   sh dev/defects-gated.sh <core|sample> <tag>
cd "$(dirname "$0")/.." || exit 1
which=$1; tag=$2
case "$which" in
  core) pkg=frmtmb; dir=tests/testthat ;;
  sample) pkg=frmtmb.sample; dir=extensions/frmtmb.sample/tests/testthat ;;
esac
out=dev/defects-log/gated-$tag
mkdir -p "$out"
export FRMTMB_BRMS_FIT_TESTS=true NOT_CRAN=true
export DEFECTS_ARM="${DEFECTS_ARM:-after}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for f in $(cat dev/defects-log/gated-$which-files.txt); do
  "$RS" dev/defects-runtest.R "$pkg" "$dir/$f" > "$out/${f%.R}.txt" 2>&1 &
done
wait
grep -h "^RESULT" $out/*.txt > dev/defects-log/gated-$tag.txt
echo "GATED $pkg ran $(wc -l < dev/defects-log/gated-$tag.txt) of $(wc -l < dev/defects-log/gated-$which-files.txt) files" >> dev/defects-log/gated-$tag.txt
