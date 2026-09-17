#!/usr/bin/env bash
# Run test files one per R process, P at a time, against a library arm.
#   dev/brmsnames-runset.sh <lane|base> <P> <outdir> file...
arm=$1; P=$2; out=$3; shift 3
mkdir -p "$out"
if [ "$arm" = lane ]; then
  LIBS="C:/Users/adf44/source/r/brmsnames-lib C:/Users/adf44/source/r/rellib-r3"
else
  LIBS="C:/Users/adf44/source/r/rellib-r3"
fi
printf '%s\n' "$@" | xargs -P "$P" -I{} bash -c '
  f="{}"; log="'"$out"'/$(echo "$f" | tr "/" "_").log"
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/brmsnames-runtest.R '"${LIBS%% *}"' "$f" '"${LIBS#* }"' > "$log" 2>&1
  echo "$(grep -h "^PASS\|^FAIL\|^ERROR\|^SKIP\|^BLOCKS" "$log" | tr -s " " | tr "\n" " ") $f"
'
