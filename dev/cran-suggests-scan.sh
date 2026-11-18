#!/usr/bin/env bash
# Lane cran: list each test file, vignette and Rd example that names a
# Suggests package (pkg::, library(pkg), data(package = "pkg")) with no
# guard for that package in the same file. A guard is
# skip_if_not_installed("pkg"), requireNamespace("pkg"), or for brms and
# rstan the skip helpers of helper-brms*.R, helper-stan.R, helper-rl.R.
# The scan is per file, not per call, so a guard anywhere in the file
# counts: it can only UNDER-report.
#
#   bash dev/cran-suggests-scan.sh > dev/cran-log/suggests-scan.txt
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
pk=$(awk '/^Suggests:/{s=1;next} /^[A-Za-z]/{s=0} s' DESCRIPTION |
  sed 's/[(,].*//; s/^ *//; s/ *$//' | grep -v '^$')
helpers='skip_unless_brms|skip_unless_brms_fit|skip_unless_brms_suite|skip_unless_rl|skip_unless_rl_stan|skip_unless_stan_identity|skip_unless_bcm'
scan() {
  local kind="$1"; shift
  for f in "$@"; do
    for p in $pk; do
      case "$p" in testthat|knitr|rmarkdown) continue;; esac
      grep -qE "\b$p::|library\($p\)|package = \"$p\"" "$f" || continue
      if grep -qE "skip_if_not_installed\(\"$p\"|requireNamespace\(\"$p\"" "$f"; then
        continue
      fi
      if [ "$p" = brms ] || [ "$p" = rstan ]; then
        grep -qE "$helpers" "$f" && continue
      fi
      n=$(grep -cE "\b$p::|library\($p\)|package = \"$p\"" "$f")
      echo "$kind $(basename "$f") $p uses=$n"
    done
  done
}
scan test tests/testthat/test-*.R tests/testthat/helper-*.R
scan vignette vignettes/*.Rmd
# examples only: the text after \examples{ in each Rd file
tmp=$(mktemp -d)
for f in man/*.Rd; do
  awk '/^\\examples[{]/{e=1} e' "$f" > "$tmp/$(basename "$f")"
  [ -s "$tmp/$(basename "$f")" ] || rm -f "$tmp/$(basename "$f")"
done
scan example "$tmp"/*.Rd
rm -rf "$tmp"
