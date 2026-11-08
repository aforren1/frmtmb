#!/usr/bin/env bash
# Review of lane loglik, attack 10: read-only scan of the other lanes'
# worktrees for (a) hunks in the files this lane edits and (b) new reads
# of logLik() or of the objective, which the split changes the meaning of.
set -u
R=C:/Users/adf44/source/r
export GIT_OPTIONAL_LOCKS=0
files="R/fit.R R/methods-fit.R R/allfit.R R/multiple.R R/par-template.R R/sugar.R R/fit-end.R R/confint.R R/refit-data.R R/compat.R R/frmtmb-package.R NAMESPACE"
echo "== loglik hunks"
for f in $files; do
  git -C $R/frmtmb-wt-loglik diff -U0 -- "$f" 2>/dev/null | grep "^@@" | sed "s|^|$f |"
done
for l in nestquad bigx nlstart superassign arstep; do
  echo "== $l hunks in files loglik edits"
  for f in $files; do
    git -C $R/frmtmb-wt-$l diff -U0 -- "$f" 2>/dev/null | grep "^@@" | sed "s|^|$f |"
  done
  echo "== $l new lines reading logLik() or the objective (R/ only)"
  git -C $R/frmtmb-wt-$l diff -- R/ 2>/dev/null | grep -n "^+.*\(logLik(\|objective\)" | grep -v "^[0-9]*:+#'"
  for u in $(git -C $R/frmtmb-wt-$l status --short -- R/ 2>/dev/null | awk '/^\?\?/ {print $2}'); do
    echo "   untracked $u:"
    grep -n "logLik(\|neg_log_prior" "$R/frmtmb-wt-$l/$u" | grep -v "#'" | head -20
  done
done
echo "== nlstart nl_sign_probe, the reads that decide"
awk '/^nl_sign_probe <- function/,/^}/' $R/frmtmb-wt-nlstart/R/par-template.R | grep -n "objective\|logLik\|prior"
echo done
