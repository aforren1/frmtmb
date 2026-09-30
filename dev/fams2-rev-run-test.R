# Reviewer's runner: one test file, one arm, one R process.
#   Rscript dev/fams2-rev-run-test.R <arm> <pkg> <file>
# <arm> is "base" (rellib-r3 with the 1f40800d test sources) or "lane"
# (the worker's library with the worktree's test sources).
args <- commandArgs(TRUE)
arm <- args[1]
pkg <- args[2]
file <- args[3]
sp <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
             "c--Users-adf44-source-r-frmtmb/",
             "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad")
root <- "C:/Users/adf44/source/r/frmtmb-wt-fams2"
base_lib <- c("C:/Users/adf44/source/r/rellib-r3",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") {
  .libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", base_lib))
  tdir <- if (pkg == "frmtmb") file.path(root, "tests/testthat") else
    file.path(root, "extensions", pkg, "tests/testthat")
} else if (arm == "r1") {
  # the round-1 state rebuilt by the reviewer (families.R as reviewed)
  # with the lane's current test sources
  .libPaths(c("C:/Users/adf44/source/r/wt-fams2-rev-lib", base_lib))
  tdir <- if (pkg == "frmtmb") file.path(root, "tests/testthat") else
    file.path(root, "extensions", pkg, "tests/testthat")
} else if (arm == "baselib") {
  # the lane's test sources against the base build: the seen-to-fail arm
  .libPaths(base_lib)
  tdir <- if (pkg == "frmtmb") file.path(root, "tests/testthat") else
    file.path(root, "extensions", pkg, "tests/testthat")
} else {
  .libPaths(base_lib)
  tdir <- file.path(sp, "base", if (pkg == "frmtmb") "core" else "sample")
}
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(sp, "stan-cache"))
options(frmtmb.brms_lp_report = TRUE)
suppressPackageStartupMessages(library(testthat))
cat("arm:", arm, " lib:", find.package(pkg), " tdir:", tdir, "\n")
res <- testthat::test_file(file.path(tdir, file), package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "summary",
                           load_package = "installed")
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s %s pass=%d fail=%d err=%d skip=%d warn=%d blocks=%d\n",
            arm, pkg, file, sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning), nrow(df)))
