# Which warning escapes test-nl-family.R on the merge, and on nlfix's own build?
#   Rscript dev/review-predread-r1-warn.R <core lib> <test tree>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/rev-predread-src/stan-cache",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
setwd(a[2])
r <- test_file("test-nl-family.R", package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(r)
cat("lib", dirname(find.package("frmtmb")), " pass", sum(d$passed), " fail", sum(d$failed), " warn", sum(d$warning), "\n")
for (t in r) for (x in t$results) if (inherits(x, "expectation_warning")) cat("WARN in", t$test, ":", substr(conditionMessage(x), 1, 200), "\n")
