#!/usr/bin/env bash
# Lane cran: R CMD build, then R CMD check --as-cran the way CRAN runs
# it, with every output line stamped by wall-clock seconds so each step's
# time can be read even below R's own [Ns] threshold.
#
#   bash dev/cran-check.sh
#
# Why NOT_CRAN is unset: CRAN never sets it, so skip_on_cran() applies.
# Why R_LIBS holds only the lane library and the user library: CRAN has
# no frmtmb extension installed, and core's CI installs none either. The
# lane library carries frmtmb (from this worktree) and drmTMB (a CRAN
# package, copied from rellib-r22, which the user library lacks).
# OMP_THREAD_LIMIT and the CPU-to-elapsed threshold copy CRAN's
# incoming Debian check (R Packages, 2e, "Releasing to CRAN"; CRAN
# policy "Checking packages ... at most 2 cores").
# _R_CHECK_TIMINGS_=0 makes R print the elapsed time of every step;
# --as-cran sets it to 10 only when it is unset. On Windows R prints
# elapsed time only, never CPU.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/dev/cran-check"
LOGD="$ROOT/dev/cran-log"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
unset NOT_CRAN FRMTMB_BRMS_FIT_TESTS FRMTMB_STAN_CACHE FRMTMB_SCALE_TESTS
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
RT=/c/rtools45
export PATH="$PANDOC:$TEX:$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC='C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools'
export R_LIBS="C:/Users/adf44/source/r/wt-cran-lib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export _R_CHECK_TIMINGS_=0
export OMP_THREAD_LIMIT=2
export _R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5
export _R_CHECK_VIGNETTE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5
stamp() { while IFS= read -r l; do printf '%s %s\n' "$(date +%s)" "$l"; done; }
rm -rf "$OUT"; mkdir -p "$OUT" "$LOGD"
cd "$OUT" || exit 1
echo "$(date +%s) BUILD START" > "$LOGD/build.log"
"$R" CMD build "$ROOT" 2>&1 | stamp >> "$LOGD/build.log"
echo "$(date +%s) BUILD END" >> "$LOGD/build.log"
tgz=$(ls -t "$OUT"/frmtmb_*.tar.gz | head -1)
[ -f "$tgz" ] || { echo "NO TARBALL" >> "$LOGD/build.log"; exit 1; }
ls -l "$tgz" >> "$LOGD/build.log"
echo "$(date +%s) CHECK START" > "$LOGD/check.log"
"$R" CMD check --as-cran --timings "$tgz" 2>&1 | stamp >> "$LOGD/check.log"
echo "$(date +%s) CHECK END" >> "$LOGD/check.log"
RC="$OUT/frmtmb.Rcheck"
cp "$RC/00check.log" "$LOGD/00check.log"
cp "$RC/00install.out" "$LOGD/00install.out"
cp "$RC/frmtmb-Ex.timings" "$LOGD/frmtmb-Ex.timings" 2>/dev/null
cp "$RC/frmtmb-Ex.Rout" "$LOGD/frmtmb-Ex.Rout" 2>/dev/null
cp "$RC/tests/testthat.Rout" "$LOGD/testthat.Rout" 2>/dev/null
cp "$RC/tests/testthat.Rout.fail" "$LOGD/testthat.Rout.fail" 2>/dev/null
# installed size by directory, in KiB, for the installed-size NOTE
( cd "$RC/frmtmb" && du -sk . * | sort -rn ) > "$LOGD/installed-du.txt"
( cd "$RC/frmtmb" && du -sk doc/* help/* R/* 2>/dev/null | sort -rn | head -40 ) \
  >> "$LOGD/installed-du.txt"
echo "CRAN-CHECK DONE" >> "$LOGD/check.log"
