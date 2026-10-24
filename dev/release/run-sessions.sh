#!/usr/bin/env bash
# Run each package's whole suite in ONE R session, three ways: in
# alphabetical order, in reverse order, and with every extension loaded
# before the first file (preload). N sessions run at a time, through
# dev/release/run-session.R, and the driver assembles one log.
#
#   bash dev/release/run-sessions.sh <tier> [N] [package ...]
#
# With no package named, all nine run, so 27 sessions. TEST_ROOT, when
# set, is the tree whose tests run (for a "before" arm read from another
# checkout); the default is this one. The caller sets the tier's
# environment (NOT_CRAN, gates), as for run-par.sh, and FRMTMB_RELLIB
# for a lane's build.
# Output: dev/release/<tier>-sessions/<package>--<order>.txt per session
# and dev/release/<tier>.log with each session's SESSION line, its
# failing DETAIL lines and its LEAK lines, then "<TIER> ran X of Y".
#
# Why beside run-par.sh: that driver isolates every file in its own
# process, which is the right instrument for a file's own result and
# the wrong one for what a file leaves behind. covr and R CMD check run
# one session per suite, and 0.71.0's coverage job failed there and
# nowhere else. Run both at a release.
set -u
tier="$1"; N="${2:-8}"; shift; [ $# -gt 0 ] && shift
# A package name in the N slot would test false in every comparison
# below, and every session would start at once
case $N in ""|*[!0-9]*) echo "N must be a number, got: $N"; exit 1;; esac
ORDERS="forward reverse preload"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TR="${TEST_ROOT:-$ROOT}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/release/run-session.R"
OUT="$ROOT/dev/release/$tier-sessions"
LOG="$ROOT/dev/release/$tier.log"
[ -f "$RUN" ] || { echo "runner missing: $RUN"; exit 1; }
pk=("$@")
[ ${#pk[@]} -gt 0 ] || pk=(frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent
                           frmtmb.learn frmtmb.ode frmtmb.sample
                           frmtmb.survival
                           frmtmb.spatial)
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
MK=C:/Users/adf44/Documents/.R/Makevars.win
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-$MK}"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
for p in "${pk[@]}"; do
  if [ "$p" = frmtmb ]; then d="$TR/tests/testthat"
  else d="$TR/extensions/$p/tests/testthat"; fi
  [ -d "$d" ] || { echo "test dir missing: $d"; exit 1; }
  # core first: it is the longest session, so it starts first
  for o in $ORDERS; do
    while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
    ( "$RS" "$RUN" "$p" "$o" "$d" > "$OUT/$p--$o.txt" 2>&1 ) &
  done
done
wait
# The log is written in one redirect, once every session has exited: a
# log appended line by line came out with two lines cut mid-line once
# (lane testiso, 2026-10-14), and the cause was not found
want=0; ran=0
up=$(echo "$tier" | tr '[:lower:]' '[:upper:]')
{
  for p in "${pk[@]}"; do
    for o in $ORDERS; do
      want=$((want + 1))
      f="$OUT/$p--$o.txt"
      hit=$(grep -a "^SESSION" "$f" | tail -1)
      if [ -n "$hit" ]; then
        ran=$((ran + 1)); echo "$hit"
        grep -a -E "^DETAIL (FAIL|ERROR|WARN) |^LEAK " "$f"
      else
        echo "NO SESSION LINE for $p $o"; tail -5 "$f"
      fi
    done
  done
  echo "$up ran $ran of $want"
} > "$LOG"
tail -1 "$LOG"
