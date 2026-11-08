#!/bin/sh
# usage: review-bigx-runjobs.sh <jobs file> <parallel>
cd "C:/Users/adf44/source/r/frmtmb-wt-bigx"
port=32000
while read lib tf name; do
  port=$((port+1))
  echo "$lib $tf $name $port"
done < "$1" | xargs -P "$2" -L 1 sh -c 'R_PARALLEL_PORT=$3 "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/review-bigx-runtest.R "$0" "$1" > "dev/review-bigx-t-$2.log" 2>&1; echo "$2 done"'
