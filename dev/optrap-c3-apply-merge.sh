#!/usr/bin/env bash
# Lane optrap, consolidation fixes: rebuild the scratch merge's R/fit.R,
# R/autoscale.R and tests/testthat/test-optrap.R from the release tree
# (frmtmb-wt-release, read only) plus this lane's diffs since build 27
# (dev/optrap-c3-base27-fit.R and dev/optrap-log/c3/base27-autoscale.R
# are build 27's sources), then install the merge.
#   bash dev/optrap-c3-apply-merge.sh
set -eu
cd "$(dirname "$0")/.."
REL=C:/Users/adf44/source/r/frmtmb-wt-release
M=dev/optrap-relmerge
T="${TEMP:-/c/Users/adf44/AppData/Local/Temp}"
for f in fit autoscale; do
  base=dev/optrap-log/c3/base27-$f.R
  tr -d '\r' < "$base" > "$T/c3-b27-$f.R"
  tr -d '\r' < "R/$f.R" > "$T/c3-new-$f.R"
  diff -u "$T/c3-b27-$f.R" "$T/c3-new-$f.R" > "dev/optrap-log/c3/$f.patch" || true
  tr -d '\r' < "$REL/R/$f.R" > "$T/c3-rel-$f.R"
  patch "$T/c3-rel-$f.R" "dev/optrap-log/c3/$f.patch"
  cp "$T/c3-rel-$f.R" "$M/R/$f.R"
done
cp tests/testthat/test-optrap.R "$M/tests/testthat/test-optrap.R"
bash dev/optrap-c3-install-merge.sh > dev/optrap-log/c3/install-merge.txt 2>&1
tail -n 1 dev/optrap-log/c3/install-merge.txt
