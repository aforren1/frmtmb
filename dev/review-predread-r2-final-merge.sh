#!/usr/bin/env bash
# Short check of fix round 2: the resolve again on latres and nlfix as finished, then the
# release fixes F1 to F4 on the merged source, each edit matching once, installed in
# rev-predread-lib/merge-final, with the files either touches.
set -u
P=/c/Users/adf44/source/r/frmtmb-wt-predread
cd "$P" || exit 1
O=dev/review-predread-out/r2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
git -C ../frmtmb-wt-nlfix diff 2>/dev/null | md5sum | cut -c1-32 > $O/merge-final-nlfix-md5.txt
bash dev/review-predread-merge2.sh > $O/merge-final-merge.log 2>&1
"$RS" dev/review-predread-merge2-resolve.R > $O/merge-final-resolve.log 2>&1
"$RS" dev/review-predread-r2-final-merge.R > $O/merge-final-edit.log 2>&1
echo finished > $O/merge-final-build.done
