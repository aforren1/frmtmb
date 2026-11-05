#!/usr/bin/env bash
# Fills the findings' generated parts from the logs: the equiv block and
# section 9. Run from the worktree root, after the runs.
f=dev/refuse-findings.md
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/refuse-gen-tables.R > /dev/null 2>&1
{
  echo 'FRMTMB_EQUIV_LIB with absolute paths; dev/refuse-f2-final.sh'
  echo 'dev/refuse-equiv-out, all 172 pairs, against the baseline'
  echo 'frmtmb-wt-r16base/dev/r16base-equiv-out (dev/refuse-verdict-diff.R,'
  echo 'dev/refuse-log/verdict-diff.txt):'
  echo
  echo '```'
  grep -a '^p[0-9]' dev/refuse-log/verdict-diff.txt | cut -c1-110
  grep -a 'headline' dev/refuse-log/verdict-diff.txt | cut -c1-90
  echo '```'
  echo
  echo 'The other 167 keep verdict and level.'
  echo
  echo '- p043 and p169 stay `defect` only because'
  echo '  `dev/equiv-data/adjudication.tsv` (rows 17 and 8) sets the verdict by'
  echo '  hand. Their per-pair files read `equivalent` and `compat row wrong` at L1.'
  echo '  **Drop both rows at consolidation.** The review ran that: 85.5.'
  echo '- p169 reads `compat row wrong`, and p010, p011, p059, p146 stay there,'
  echo '  because `sample.tsv` freezes the status the flow compares against; the'
  echo '  rows themselves are fixed.'
  echo '- p029 moved because the translator reads the univariate `Ymi`. It and p169'
  echo '  are the complete-case route: brms on the observed rows.'
} > dev/refuse-log/equiv-block.md
awk 'FNR==NR { b = b $0 "\n"; next } /^REFUSE-EQUIV-BLOCK$/ { printf "%s", b; next } { print }' \
  dev/refuse-log/equiv-block.md "$f" > dev/refuse-log/findings.tmp
sed '/^## 9\. Generated tables and run counts/,$d' dev/refuse-log/findings.tmp > "$f"
{
  echo "## 9. Generated tables and run counts"
  echo
  cat dev/refuse-log/gen-tables.md
  echo
  echo "### mi() predictor under each family (dev/refuse-mi-pred.R, seed 170)"
  echo
  echo '```'
  grep -aE 'frm:' dev/refuse-log/mi-pred.log | cut -c1-90
  echo '```'
  echo
  echo "### B1 and the read arm, per family (dev/refuse-f1-unread.R)"
  echo
  echo '```'
  grep -ahE '^(family|UNREAD|READ)' dev/refuse-f1-out/{gaussian,student,skew_normal,beta,Gamma,weibull,exponential,frechet,inverse.gaussian,lognormal,exgaussian,gen_extreme_value,von_mises,asym_laplace}.txt | cut -c1-110
  echo '```'
  echo
  echo "### Tests on the final build"
  echo
  echo '```'
  awk '/./' dev/refuse-log/final-counts.txt
  echo '```'
} >> "$f"
