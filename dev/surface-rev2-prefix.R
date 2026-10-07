# Reviewer re-check: run a punch-round test file against the PRE-FIX
# build of the first review (scratchpad installs of the sources the
# first review installed), to see the new assertions fail.
#   Rscript dev/surface-rev2-prefix.R oldsample|oldcore <pkg> <file>
a <- commandArgs(TRUE)
S <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
            "c--Users-adf44-source-r-frmtmb/",
            "7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/rev2-", a[1])
.libPaths(c(S, "C:/Users/adf44/source/r/surface-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
pkg <- a[2]
suppressPackageStartupMessages({
  library(testthat); library(pkg, character.only = TRUE)
})
cat("frmtmb", find.package("frmtmb"), "| frmtmb.sample",
    tryCatch(find.package("frmtmb.sample"), error = function(e) "-"), "\n")
tdir <- if (pkg == "frmtmb") "tests/testthat" else
  file.path("extensions", pkg, "tests/testthat")
res <- test_file(file.path(tdir, a[3]), package = pkg, env = test_env(pkg),
                 reporter = "silent", stop_on_failure = FALSE)
for (t in res) for (x in t$results) {
  if (inherits(x, c("expectation_failure", "expectation_error"))) {
    cat("FAIL [", t$test, "]:", substr(gsub("\n", " ", conditionMessage(x)),
                                      1, 220), "\n")
  }
}
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d\n", a[1], a[3],
            sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped)))
