#!/usr/bin/env bash
# After the restart of 2026-11-02: reinstall, rerun what R/compat.R's last note
# edit reaches and what the restart cut off (three R CMD checks).
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
LIB=C:/Users/adf44/source/r/wt-nlfix-lib
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=true
step() { echo "$(date +%H:%M:%S) $*" >> dev/nlfix-final5.progress; }
: > dev/nlfix-final5.progress
step install
bash dev/nlfix-install.sh > dev/nlfix-install-final.log 2>&1
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for p in survival eam learn; do
  "$RC" CMD INSTALL --library="$LIB" --no-multiarch extensions/frmtmb.$p > dev/nlfix-install-$p.log 2>&1 || echo "INSTALL FAILED $p" >> dev/nlfix-final5.progress
done
unset R_LIBS
step "tests the note reaches"
rm -f dev/nlfix-suite-core/frmtmb--test-compat*.txt dev/nlfix-suite-core/frmtmb--test-nl*.txt
ONLY='test-compat|test-nl' bash dev/nlfix-suite.sh core 0 8 frmtmb > dev/nlfix-suite-core.log 2>&1
step "R CMD check"
bash dev/nlfix-check.sh . > dev/nlfix-check-frmtmb.out 2>&1 &
bash dev/nlfix-check.sh extensions/frmtmb.sample > dev/nlfix-check-frmtmb.sample.out 2>&1 &
bash dev/nlfix-check.sh extensions/frmtmb.eam > dev/nlfix-check-frmtmb.eam.out 2>&1 &
wait
step done
echo FINAL5 DONE > dev/nlfix-final5.done
