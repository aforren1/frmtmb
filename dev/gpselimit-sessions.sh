#!/usr/bin/env bash
# Lane gpselimit: whole suites in ONE R session each, through
# dev/release/run-session.R, with the lane's library first. Each argument
# is "<package>:<order>"; two sessions at a time. Logs go to
# dev/gpselimit-out/sessions/<package>--<order>.txt, and the SESSION,
# failing DETAIL and LEAK lines to dev/gpselimit-out/sessions.log.
# Usage: bash dev/gpselimit-sessions.sh frmtmb:forward frmtmb.sample:reverse
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-gpselimit
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=/c/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-out/sessions
LOG=/c/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-out/sessions.log
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-gpselimit-lib
mkdir -p "$OUT"
for job in "$@"; do
  p=${job%%:*}; o=${job##*:}
  if [ "$p" = frmtmb ]; then d="$WT/tests/testthat"
  else d="$WT/extensions/$p/tests/testthat"; fi
  while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 1; done
  ( "$RS" "$WT/dev/release/run-session.R" "$p" "$o" "$d" \
      > "$OUT/$p--$o.txt" 2>&1 ) &
done
wait
{
  for job in "$@"; do
    p=${job%%:*}; o=${job##*:}
    f="$OUT/$p--$o.txt"
    grep -a "^lib" "$f" | head -1
    hit=$(grep -a "^SESSION" "$f" | tail -1)
    if [ -n "$hit" ]; then echo "$p $o $hit"
    else echo "NO SESSION LINE for $p $o"; tail -5 "$f"; fi
    grep -a -E "^DETAIL (FAIL|ERROR|WARN) |^LEAK " "$f"
  done
} >> "$LOG"
