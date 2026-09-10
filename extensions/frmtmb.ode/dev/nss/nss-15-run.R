Sys.setenv(NOT_CRAN = "true")
library(testthat)
library(frmtmb)
library(frmtmb.ode)
cat("frmtmb.ode from", dirname(find.package("frmtmb.ode")), "\n")
setwd(paste0("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/",
             "frmtmb.ode/tests/testthat"))
r <- test_file("test-ode-events.R", reporter = "summary",
               package = "frmtmb.ode")
d <- as.data.frame(r)
cat("\nFILE test-ode-events.R files=1 pass=", sum(d$passed),
    " fail=", sum(d$failed), " err=", sum(d$error),
    " skip=", sum(d$skipped), " nb_tests=", nrow(d), "\n", sep = "")
bad <- d[d$failed > 0 | d$error > 0, "test"]
if (length(bad)) cat("BROKE:", paste(bad, collapse = " | "), "\n")
