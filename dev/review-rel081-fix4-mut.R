# Release review of 0.81.0, confirmation of fix4: does
# test-loo-newdata-fitted.R see the pp_check() check go? The
# loo_check_finite() line of draws_ppc_psis() is removed in the loaded
# namespace, then the file runs.
#   Rscript dev/review-rel081-fix4-mut.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb.sample")), "sample",
    format(packageVersion("frmtmb.sample")), "\n")
ns <- asNamespace("frmtmb.sample")
f <- get("draws_ppc_psis", ns)
src <- deparse(f)
hit <- grep("loo_check_finite(", src, fixed = TRUE)
stopifnot(length(hit) == 1L)
src <- src[-hit]
g <- eval(parse(text = src))
environment(g) <- environment(f)
unlockBinding("draws_ppc_psis", ns)
assign("draws_ppc_psis", g, envir = ns)
lockBinding("draws_ppc_psis", ns)
cat("draws_ppc_psis() without loo_check_finite()\n")
tf <- "C:/Users/adf44/source/r/frmtmb-wt-rel081/extensions/frmtmb.sample/tests/testthat/test-loo-newdata-fitted.R"
res <- as.data.frame(test_file(tf, package = "frmtmb.sample",
                               env = test_env("frmtmb.sample"),
                               reporter = "silent", stop_on_failure = FALSE))
print(res[c("test", "nb", "failed", "error", "skipped")], row.names = FALSE)
cat("total: tests", nrow(res), "failed", sum(res$failed), "errors",
    sum(res$error), "\n")
cat("FIX4 MUT DONE\n")
