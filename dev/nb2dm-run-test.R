# One test file, one R process, against a chosen library.
#   Rscript dev/nb2dm-run-test.R <package> <test file> [lib]
# lib defaults to the lane's private library; "base" means rellib-r10.
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
lib <- if (length(a) >= 3L) a[3] else "C:/Users/adf44/source/r/wt-nb2dm-lib"
base <- "C:/Users/adf44/source/r/rellib-r10"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(unique(if (identical(lib, "base")) c(base, user) else
  c(lib, base, user)))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) Sys.setenv(R_MAKEVARS_USER = mk)
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(utils::packageVersion("frmtmb")), p,
    dirname(system.file(package = p)), format(utils::packageVersion(p)),
    "\n")
r <- as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " warn=", sum(r$warning), " blocks=", nrow(r), "\n", sep = "")
fails <- r[r$failed > 0 | r$error, "test"]
if (length(fails)) cat("FAILING BLOCK:", fails, sep = "\n  ")
