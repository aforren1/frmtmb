# Lane sampfront: run ONE test file of a package against a library stack.
#   Rscript dev/sampfront-run-test.R <pkg> <test file> [base]
# With "base" the lane library is dropped (the "before" arm, rellib-r20).
local({
args <- commandArgs(TRUE)
pkg <- args[1]
f <- args[2]
base <- length(args) > 2 && identical(args[3], "base")
lib <- "C:/Users/adf44/source/r/wt-sampfront-lib"
.libPaths(c(if (!base) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampfront/dev/stan-cache")
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib: frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
if (pkg != "frmtmb") {
  cat("lib:", pkg, format(packageVersion(pkg)), find.package(pkg), "\n")
}
cat("gated:", Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "\n")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
for (r in res) {
  for (e in r$results) {
    k <- if (inherits(e, "expectation_failure")) "FAIL" else
      if (inherits(e, "expectation_error")) "ERROR" else
        if (inherits(e, "expectation_skip")) "SKIP" else
          if (inherits(e, "expectation_warning")) "WARN" else NA
    if (!is.na(k)) {
      cat("DETAIL ", k, " [", r$test, "] ",
          substr(gsub("\n", " | ", conditionMessage(e)), 1, 700), "\n",
          sep = "")
    }
  }
}
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d %.1fs\n",
            pkg, basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning),
            proc.time()[["elapsed"]] - t0))
})
