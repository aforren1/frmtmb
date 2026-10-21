# One test file, one R process, against the lane library or the base.
#   Rscript dev/loodraws-runtest.R <package> <test file> [lane|base]
# Prints the lib: line and run-tests.R's RESULT line, then every
# failing, erroring or skipped block and every escaped warning.
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
arm <- if (length(a) >= 3L) a[3] else "lane"
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(utils::packageVersion("frmtmb")), "|",
    dirname(system.file(package = p)), p,
    format(utils::packageVersion(p)), "\n")
r <- as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                             reporter = "check", stop_on_failure = FALSE))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " warn=", sum(r$warning), " blocks=", nrow(r), "\n", sep = "")
bad <- r[r$failed > 0 | r$error | r$skipped | r$warning > 0, "test"]
if (length(bad)) cat("NOT CLEAN:", bad, sep = "\n  ")
