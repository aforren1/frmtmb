#!/bin/sh
# Reviewer: install the lane's core and frmtmb.sample, unchanged, into the
# reviewer's own library; then roxygenise a COPY of the tree (never the
# worktree) to show whether man/ and NAMESPACE are current.
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-nb2dm
LIB=C:/Users/adf44/source/r/wt-nb2dm-rev-lib
OUT=$WT/dev/nb2dm-rev-out
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"$RC" CMD INSTALL --library=$LIB --no-multiarch $WT > $OUT/install-core.txt 2>&1
"$RC" CMD INSTALL --library=$LIB --no-multiarch $WT/extensions/frmtmb.sample > $OUT/install-sample.txt 2>&1
CP=C:/Users/adf44/AppData/Local/Temp/3/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/nb2dm-rev-rox
rm -rf $CP; mkdir -p $CP/core $CP/sample
cp -r $WT/R $WT/man $WT/NAMESPACE $WT/DESCRIPTION $CP/core/
cp -r $WT/extensions/frmtmb.sample/R $WT/extensions/frmtmb.sample/man $WT/extensions/frmtmb.sample/NAMESPACE $WT/extensions/frmtmb.sample/DESCRIPTION $CP/sample/
for p in core sample; do
  "$R" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r10','C:/Users/adf44/AppData/Local/R/win-library/4.6')); suppressMessages(roxygen2::roxygenise('$CP/$p'))" > $OUT/rox-$p.txt 2>&1 || true
done
diff -r $WT/man $CP/core/man > $OUT/rox-diff.txt 2>&1 || true
diff $WT/NAMESPACE $CP/core/NAMESPACE >> $OUT/rox-diff.txt 2>&1 || true
diff -r $WT/extensions/frmtmb.sample/man $CP/sample/man >> $OUT/rox-diff.txt 2>&1 || true
echo "rox diff lines: $(wc -l < $OUT/rox-diff.txt)" >> $OUT/rox-diff.txt
