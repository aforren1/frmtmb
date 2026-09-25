#!/bin/bash
# Reviewer 2: one worker of the item-6 coverage run. Waits for 5 GB free
# before each seed. Usage: r2-cov-driver.sh <lib> <tag> <seed...>
cd /c/Users/adf44/source/r/frmtmb-wt-phase3b
lib=$1; tag=$2; shift 2
log=dev/phase3b-review2/cov-$lib-$tag.log
for s in "$@"; do
  [ -f dev/phase3b-review2/cov-$lib/seed-$s.rds ] && continue
  while :; do
    fr=$(powershell -NoProfile -Command "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB)")
    fr=$(echo "$fr" | tr -d '\r')
    [ "$fr" -ge 5 ] && break
    echo "$(date +%T) waiting, free ${fr} GB" >> $log
    sleep 60
  done
  echo "$(date +%T) start $s free ${fr} GB" >> $log
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/phase3b-review2/r2-cov-base.R $lib $s >> $log 2>&1
done
echo "$(date +%T) DONE" >> $log
