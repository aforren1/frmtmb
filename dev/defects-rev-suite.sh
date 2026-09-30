#!/bin/sh
# Reviewer of lane defects: every test file of a package (or the files
# named after the tag), one R process each, P at a time.
#   P=16 sh dev/defects-rev-suite.sh <arm> <pkg> <tag> [file ...]
# Output: dev/defects-rev-log/<tag>/<file>.txt and <tag>.txt (RESULT lines).
W=C:/Users/adf44/source/r/frmtmb-wt-defects
cd "$W" || exit 1
arm=$1; pkg=$2; tag=$3; shift 3; P=${P:-16}
case "$pkg" in
  frmtmb) dir=tests/testthat ;;
  *) dir=extensions/$pkg/tests/testthat ;;
esac
if [ $# -gt 0 ]; then files="$*"; else files=$(cd $dir; ls test-*.R); fi
out=dev/defects-rev-log/$tag
mkdir -p "$out"
export FRMTMB_STAN_CACHE=C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/defrev-stan-cache
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for k in $(seq 1 "$P"); do
  (
    j=0
    for f in $files; do
      j=$((j + 1))
      [ $(( (j - 1) % P + 1 )) -eq "$k" ] || continue
      b=$(basename "$f" .R)
      (cd $dir && "$RS" "$W/dev/defects-rev-runtest.R" "$arm" "$pkg" "$f" \
        > "$W/$out/$b.txt" 2>&1)
    done
  ) &
done
wait
n=$(echo $files | wc -w)
grep -h "^RESULT" $out/*.txt > dev/defects-rev-log/$tag.txt
got=$(wc -l < dev/defects-rev-log/$tag.txt)
echo "SUITE $arm $pkg ran $got of $n files" >> dev/defects-rev-log/$tag.txt
for f in $out/*.txt; do
  grep -q "^RESULT" "$f" || echo "NO RESULT $(basename "$f")" >> dev/defects-rev-log/$tag.txt
done
