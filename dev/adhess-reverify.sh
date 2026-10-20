#!/usr/bin/env bash
# Lane adhess, resumed 2026-10-19: rerun the measurements of sections 1
# to 6 and the new test file on the fresh install, lane and base. Logs
# in dev/adhess-log/rv/. At most N R processes at once.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-adhess
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS
L=dev/adhess-log/rv
mkdir -p $L
N=${N:-8}
jw() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
go() { # <log> <arm> <script> args...
  local lg=$1 arm=$2; shift 2
  jw
  echo "$(date +%T) START $lg" >> $L/driver.txt
  ( ADHESS_ARM=$arm "$RS" "$@" < /dev/null > $L/$lg.txt 2>&1
    echo "$(date +%T) END $lg rc $?" >> $L/driver.txt ) &
}
LIB=C:/Users/adf44/source/r/wt-adhess-lib
TF=tests/testthat/test-se-exact-hessian.R
go t3-exact-lane lane dev/ciharden-run1.R $LIB frmtmb $TF
go t3-exact-base base dev/ciharden-run1.R base frmtmb $TF
go tolmargin-lane lane dev/adhess-tolmargin.R
go tolmargin-base base dev/adhess-tolmargin.R
for a in lane base; do
  go restep-$a $a dev/adhess-restep.R
  go move-$a $a dev/adhess-move.R $a $L/move-$a.rds
  go fdcov-$a $a dev/adhess-fdcov.R
  go glm-$a $a dev/adhess-glm.R
  go nls-$a $a dev/adhess-nls.R
  go shift-$a $a dev/adhess-shift.R
  go routes-$a $a dev/adhess-routes.R
done
wait
"$RS" dev/adhess-movecmp.R $L/move-lane.rds $L/move-base.rds > $L/movecmp.txt 2>&1
echo "$(date +%T) all done" >> $L/driver.txt
