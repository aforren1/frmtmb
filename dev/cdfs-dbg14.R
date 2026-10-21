# Is asym_laplace's false convergence the censoring's or the family's?
# test-cdfs.R's data, uncensored, on the base build and on this lane's.
a <- commandArgs(TRUE)
.libPaths(c(if (length(a)) a[1], "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(43); n <- 160; x <- rnorm(n)
y <- brms::rskew_normal(n, 1 + 0.5 * x, 1, 4); y <- brms::rexgaussian(n, 2 + 0.5 * x, 0.5, 1)
y <- brms::rasym_laplace(n, 1 + 0.5 * x, 1, 0.3)
d <- data.frame(y = y, x = x)
fit <- withCallingHandlers(frm(y ~ x, family = asym_laplace(), data = d), warning = function(w) { cat("WARN:", substr(conditionMessage(w), 1, 60), "\n"); invokeRestart("muffleWarning") })
cat(dirname(find.package("frmtmb")), "uncensored conv", fit$opt$convergence, fit$opt$message, "\n")
