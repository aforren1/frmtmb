#!/usr/bin/env bash
# Run test files one per R process on this lane's library.
#   bash dev/refuse-runtests.sh <outdir under dev/> <jobs> <job list file>
# job list: "<package> <test file>" per line. Prints the RESULT lines.
out="dev/$1"; N="$2"; list="$3"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export NOT_CRAN=true
mkdir -p "$out"
export out R
cat "$list" | xargs -P "$N" -L1 bash -c \
  'p="$0"; f="$1"; timeout 3000 "$R" dev/refuse-runtest.R "$p" "$f" > "$out/$p--$(basename "$f" .R).txt" 2>&1'
cat "$out"/*.txt | grep -a '^RESULT' | sort
grep -al '^DETAIL\|^Error\|Execution halted' "$out"/*.txt
