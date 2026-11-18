#!/usr/bin/env bash
# Lane cran: time every vignette (one process each) and every core test
# file (one process each) with NOT_CRAN unset.
#
#   bash dev/cran-runall.sh
#
# The number of processes at once is read from dev/cran-log/par-n before
# each launch, so it can be 1 while R CMD check holds the lane's other
# process slot and 2 after (the round allows two R processes per lane).
# Output: dev/cran-log/vig/<name>.txt, dev/cran-log/tests/<file>.txt,
# and wall-clock spans in dev/cran-log/runall-spans.tsv.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOGD="$ROOT/dev/cran-log"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
unset NOT_CRAN FRMTMB_BRMS_FIT_TESTS FRMTMB_STAN_CACHE FRMTMB_SCALE_TESTS
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
RT=/c/rtools45
export PATH="$PANDOC:$TEX:$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC='C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools'
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
# what --as-cran sets for the check's child processes, and CRAN's
# incoming thread limit
export _R_CHECK_LIMIT_CORES_=TRUE
export OMP_THREAD_LIMIT=2
mkdir -p "$LOGD/vig" "$LOGD/tests"
SPANS="$LOGD/runall-spans.tsv"
: > "$SPANS"
VSCR="$ROOT/dev/cran-check/vigwork"
slots() { n=$(cat "$LOGD/par-n" 2>/dev/null); echo "${n:-1}"; }
wait_slot() {
  while [ "$(jobs -rp | wc -l)" -ge "$(slots)" ]; do sleep 2; done
}
run1() {
  local kind="$1" name="$2" out="$3"; shift 3
  local s e
  s=$(date +%s.%N)
  "$RS" "$@" > "$out" 2>&1
  e=$(date +%s.%N)
  printf '%s\t%s\t%s\t%s\n' "$kind" "$name" "$s" "$e" >> "$SPANS"
}
for v in "$ROOT"/vignettes/*.Rmd; do
  wait_slot
  b=$(basename "$v" .Rmd)
  run1 vig "$b" "$LOGD/vig/$b.txt" "$ROOT/dev/cran-vignette.R" "$v" "$VSCR" &
done
for f in "$ROOT"/tests/testthat/test-*.R; do
  wait_slot
  b=$(basename "$f" .R)
  run1 test "$b" "$LOGD/tests/$b.txt" "$ROOT/dev/cran-testfile.R" "$f" &
done
wait
echo "RUNALL DONE $(grep -c . "$SPANS") spans" >> "$LOGD/runall-done.txt"
