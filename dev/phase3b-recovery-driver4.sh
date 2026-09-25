#!/usr/bin/env bash
# Punch round 1, after the 09:24 crash: refit the missing seeds, writing
# to recov5/. The build history, because the first two launches failed:
#   phase3b-lib7: every cleft fit ended in std::bad_alloc (the
#     distribution functions were formed on all 12,000 rows);
#   phase3b-lib8: formed on the censored rows only (ddm_on_rows());
#     cleft still peaked near 8 GB per fit and 31 of 36 ran out beside
#     four others;
#   phase3b-lib9: interval rows scored at their boundary only, the
#     lower edge's function reused. One cleft fit peaks at 3.2 GB, cens
#     at 3.0 GB and contfix at 3.7 GB (dev/phase3b-peakmem.ps1).
# All three agree to 1e-13 in logLik (dev/phase3b-log/cleft-lib7-vs-
# lib8.txt, cleft-lib8-vs-lib9.txt), and the arms without left or
# interval censoring never reach the changed code. THREE workers, each
# waiting for 5 GB of free memory before each seed; the machine went
# down with eleven R processes of this lane and other lanes' beside.
# Targets in dev/phase3b-recov5-todo.R, minus the seeds
# dev/phase3b-log/rds-check.txt lists as already fitted.
cd "$(dirname "$0")/.." || exit 1
export P3B_LIB=C:/Users/adf44/source/r/phase3b-lib9
R="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/phase3b-log/recov5-logs
freegb() {
  powershell -NoProfile -Command \
    "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB)"
}
work() {  # name, then "arm seed" pairs
  local name=$1; shift
  while [ $# -gt 1 ]; do
    until [ "$(freegb)" -ge 5 ]; do sleep 30; done
    "$R" dev/phase3b-eam-recovery7.R "$1" "$2" \
      >> "dev/phase3b-log/recov5-logs/$name.log" 2>&1
    shift 2
  done
}
"$R" dev/phase3b-recov5-todo.R || exit 1
mapfile -t jobs < dev/phase3b-log/recov5-todo.txt
for w in 0 1 2; do
  args=()
  for ((i = w; i < ${#jobs[@]}; i += 3)); do args+=(${jobs[$i]}); done
  work "w$w" "${args[@]}" &
done
wait
date -Iseconds > dev/phase3b-log/recov5.done
