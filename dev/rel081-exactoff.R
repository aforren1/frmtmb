# The new exact-reader tests of test-loo-newdata-fitted.R on rellib-r20
# with the reader's one argument turned back: draws_row_loglik() calling
# row_lpdf(exact = FALSE), the fitted objective's barrier, as frmtmb.sample
# 0.27.0 read it. Each test's failures and errors are counted, so the
# record shows which tests see the argument go.
#
#   Rscript dev/rel081-exactoff.R > dev/rel081-log/exactoff.txt 2>&1
LIB <- "C:/Users/adf44/source/r/rellib-r20"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/stan-cache")
suppressMessages({library(testthat); library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb.sample")), "sample",
    format(packageVersion("frmtmb.sample")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb.sample")
f <- get("draws_row_loglik", ns)
src <- deparse(f)
hit <- grep("exact = TRUE", src, fixed = TRUE)
stopifnot(length(hit) == 1L)
src[hit] <- sub("exact = TRUE", "exact = FALSE", src[hit], fixed = TRUE)
g <- eval(parse(text = src))
environment(g) <- environment(f)
unlockBinding("draws_row_loglik", ns)
assign("draws_row_loglik", g, envir = ns)
lockBinding("draws_row_loglik", ns)
cat("draws_row_loglik() now calls row_lpdf(exact = FALSE)\n")
tf <- "C:/Users/adf44/source/r/frmtmb-wt-rel081/extensions/frmtmb.sample/tests/testthat/test-loo-newdata-fitted.R"
res <- as.data.frame(test_file(tf, package = "frmtmb.sample",
                               env = test_env("frmtmb.sample"),
                               reporter = "silent", stop_on_failure = FALSE))
print(res[c("test", "nb", "failed", "error", "skipped")], row.names = FALSE)
cat("total: tests", nrow(res), "failed", sum(res$failed), "errors",
    sum(res$error), "\n")
cat("EXACTOFF DONE\n")
