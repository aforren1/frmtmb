# Run one test file against the lane's build, one file per process.
# Usage: Rscript dev/famuni-test1.R <package> <test file> [base]
# With "base" the lane library is dropped, so the file runs against
# rellib-r8 (the 0.70.0 build) for a "before" arm.
args <- commandArgs(TRUE)
pkg <- args[1]
f <- args[2]
base <- isTRUE(args[3] == "base")
lp <- c("C:/Users/adf44/source/r/rellib-r8",
        "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (!base) lp <- c("C:/Users/adf44/source/r/wt-famuni-lib", lp)
.libPaths(lp)
Sys.setenv(NOT_CRAN = "true")
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) Sys.setenv(R_MAKEVARS_USER = mk)
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/stan-cache")
}
suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), "| frmtmb",
    format(packageVersion("frmtmb")), "|", pkg, "from", find.package(pkg),
    "\n")
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "check", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s tests=%d pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), nrow(df), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
