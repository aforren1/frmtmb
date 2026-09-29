# Lane wt-resmooth. Verify the delivered library with a FIT, not a
# version string: a hollow install still answers packageVersion().
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("installed:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), "\n")
cat("smooth_b_idx present:",
    exists("smooth_b_idx", envir = asNamespace("frmtmb")), "\n")
set.seed(5)
n <- 300
d <- data.frame(x = runif(n), g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + rnorm(10, 0, 0.5)[d$g] * d$x + rnorm(n, 0, 0.3)
f <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)), data = d))
a <- fitted(f, re_formula = NA)[, "Estimate"]
b <- fitted(f, re_formula = NULL)[, "Estimate"]
cat("logLik:", sprintf("%.10f", as.numeric(logLik(f))), "\n")
cat("max abs difference, NA against NULL:",
    sprintf("%.3e", max(abs(a - b))), "\n")
cat("simulate(NA) identical to the conditional draws:",
    identical(simulate(f, nsim = 5, seed = 1, re_formula = NA),
              simulate(f, nsim = 5, seed = 1)), "\n")
