#!/usr/bin/env bash
# bash dev/nlfix-l1-run.sh <outdir> <spec> [<spec> ...]: L1 of each spec, one
# R process each, 6 at once. Rtools at the END of PATH (dev/equiv-run.sh).
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
out=$1; shift
mkdir -p "$out/logs"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=true
printf '%s\n' "$@" | xargs -P 6 -I{} bash -c \
  "timeout 1800 \"$R\" dev/nlfix-l1.R {} \"$out\" > \"$out/logs/{}.log\" 2>&1 || echo FAILED {} >> \"$out/logs/{}.log\""
for s in "$@"; do
  f="$out/$s.txt"
  if [ -f "$f" ]; then
    printf '%-26s ' "$s"
    grep -E '\$ (verdict|l1_const|l1_grad|l1_extra_const|l1_ok|frm_error|brms_error|l1_error|l2|reason)' "$f" |
      sed -E 's/^ *\$ //; s/ *: / /' | cut -c1-120 | tr '\n' ';'
    echo
  else
    echo "$s: NO RESULT"; tail -3 "$out/logs/$s.log" | cut -c1-200
  fi
done
