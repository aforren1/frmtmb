#!/usr/bin/env bash
# Launch the recovery arms as parallel seed blocks. Each worker writes
# its own log under dev/phase3b-log/recov-logs/. Counts are taken from
# the RDS files, never from this launcher.
cd "$(dirname "$0")/.." || exit 1
R="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/phase3b-log/recov-logs
run() {  # arm from to
  "$R" dev/phase3b-eam-recovery.R "$1" "$2" "$3" \
    > "dev/phase3b-log/recov-logs/$1-$2-$3.log" 2>&1
}
runs() {  # from to
  "$R" dev/phase3b-session-recovery.R "$1" "$2" \
    > "dev/phase3b-log/recov-logs/session-$1-$2.log" 2>&1
}
case "$1" in
  eam)
    run cens 1 50 & run cens 51 100 &
    run cont 1 34 & run cont 35 67 & run cont 68 100 &
    run collapse 1 30 & run collapse 31 60 &
    run contdef 1 20 &
    wait ;;
  session)
    runs 1 101 & runs 102 202 &
    wait ;;
esac
