#!/usr/bin/env bash
# Run one extension's test suite, ONE FILE PER PROCESS, against either
# the base build (rellib-r3) or this lane's build.
# Usage: bash dev/phase3b-suite.sh <pkg> <base|new> [newlib] [label]
# Output: dev/phase3b-log/suite-<pkg>-<label>-<tag>.txt, one RESULT line per
# file plus each file's full reporter output in the matching directory.
cd "$(dirname "$0")/.." || exit 1
pkg=$1; tag=$2; newlib=${3:-C:/Users/adf44/source/r/phase3b-lib}; label=${4:-plain}
R="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/phase3b-log/suite-$pkg-$label-$tag
mkdir -p "$out"
: > "$out.txt"
if [ "$tag" = base ]; then libs=""; else
  libs="\"$newlib\", \"C:/Users/adf44/source/r/phase3b-lib\","; fi
runner="$out/runner.R"
cat > "$runner" <<EOF
.libPaths(unique(c($libs "C:/Users/adf44/source/r/rellib-r3",
  "C:/Users/adf44/source/r/pinlib",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(NOT_CRAN = "true")
a <- commandArgs(TRUE)
suppressPackageStartupMessages({library(testthat); library(a[2], character.only = TRUE)})
cat("from", find.package(a[2]), find.package("frmtmb.eam"), "\n")
cat("env:", Sys.getenv(c("FRMTMB_BRMS_FIT_TESTS", "FRMTMB_DRMTMB_FIT_TESTS",
                         "FRMTMB_FUZZ", "FRMTMB_STAN_CACHE",
                         "R_MAKEVARS_USER", "TMP")), "\n")
cat("StanHeaders", tryCatch(as.character(packageVersion("StanHeaders")),
                            error = function(e) "absent"), "\n")
r <-as.data.frame(test_file(a[1], reporter = "progress", package = a[2],
                             load_package = "installed"))
print(r[, c("test", "nb", "failed", "skipped", "error", "passed")])
cat("RESULT", basename(a[1]), "blocks", nrow(r), "expectations", sum(r\$nb),
    "pass", sum(r\$passed), "fail", sum(r\$failed), "err", sum(r\$error),
    "skip", sum(r\$skipped), "\n")
EOF
for f in extensions/$pkg/tests/testthat/test-*.R; do
  b=$(basename "$f" .R)
  (cd "extensions/$pkg" && "$R" "../../$runner" "tests/testthat/$b.R" "$pkg") \
    > "$out/$b.txt" 2>&1
  line=$(grep "^RESULT" "$out/$b.txt" || echo "RESULT $b.R NO-RESULT")
  echo "$line" >> "$out.txt"
done
echo "DONE $(date -Iseconds)" >> "$out.txt"
