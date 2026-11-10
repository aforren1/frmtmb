# Lane hollow: run ONE test file of a package against a library stack.
#   Rscript dev/hollow-run-test.R <pkg> <test file> [base]
# With "base" the lane library is dropped (the "before" arm, rellib-r19).
args <- commandArgs(TRUE)
pkg <- args[1]
f <- args[2]
base <- length(args) > 2 && identical(args[3], "base")
lib <- "C:/Users/adf44/source/r/wt-hollow-lib"
.libPaths(c(if (!base) lib, "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-hollow/dev/stan-cache")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
if (pkg != "frmtmb") {
  cat("lib:", find.package(pkg), format(packageVersion(pkg)), "\n")
}
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
fails <- 0L
for (r in res) {
  for (e in r$results) {
    if (inherits(e, c("expectation_failure", "expectation_error"))) {
      cat("FAIL/ERROR in '", r$test, "': ",
          substr(conditionMessage(e), 1, 600), "\n", sep = "")
    }
    if (inherits(e, "expectation_warning")) {
      cat("WARN in '", r$test, "': ", substr(conditionMessage(e), 1, 300),
          "\n", sep = "")
    }
    if (inherits(e, "expectation_skip")) {
      cat("SKIP in '", r$test, "': ", substr(conditionMessage(e), 1, 200),
          "\n", sep = "")
    }
  }
}
cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d %.1fs\n",
            pkg, basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning),
            proc.time()[["elapsed"]] - t0))
