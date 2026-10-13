#!/usr/bin/env bash
# Driver for the rel070 moparity reruns: at most 6 R processes at once.
set -u
REL=C:/Users/adf44/source/r/frmtmb-wt-release
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$REL/dev/rel070-log/studies/moparity"
mkdir -p "$OUT/dev/moparity-log"
export TMP="C:/Users/adf44/AppData/Local/Temp/3" TEMP="C:/Users/adf44/AppData/Local/Temp/3"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$REL/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
unset FRMTMB_RELLIB
N=6
run() { # name, cwd, args...
  local nm="$1" wd="$2"; shift 2
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  ( cd "$wd" && echo "CMD: Rscript $*" > "$OUT/$nm.txt" &&
    "$RS" "$@" >> "$OUT/$nm.txt" 2>&1; echo "EXIT $?" >> "$OUT/$nm.txt" ) &
}
run cox-study-r8 "$OUT" "$REL/dev/rel070-moparity-cox-study.R" lane 200
run test-reparam "$REL" dev/release/run-tests.R frmtmb.sample \
  extensions/frmtmb.sample/tests/testthat/test-reparam.R
run gate-newcor-0-25 "$REL" dev/rel070-moparity-rev2-gate.R newcor 0 25 32
run gate-newcor-1-13 "$REL" dev/rel070-moparity-rev2-gate.R newcor 1 13 16
run gate-newcor-0-33 "$REL" dev/rel070-moparity-rev2-gate.R newcor 0 33 33
run gate-newcor-1-33 "$REL" dev/rel070-moparity-rev2-gate.R newcor 1 33 33
run update "$REL" dev/rel070-moparity-rev3-update.R
run rev2-update "$REL" dev/rel070-moparity-rev2-update.R
run names-r8 "$REL" dev/moparity-rev3-names.R C:/Users/adf44/source/r/rellib-r8
run nl-r8 "$REL" dev/moparity-rev3-nl.R C:/Users/adf44/source/r/rellib-r8
wait
echo DONE
