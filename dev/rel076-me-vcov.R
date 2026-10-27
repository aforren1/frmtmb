# Does marginaleffects read a user's vcov = matrix for a frmtmb fit?
# The guard's inverse case is a doubled covariance, which must double
# the variance (sqrt(2) on the standard error) if the matrix is read.
#   Rscript dev/rel076-me-vcov.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "marginaleffects",
    as.character(packageVersion("marginaleffects")), "\n")
set.seed(1)
d <- data.frame(x = stats::rnorm(80))
d$y <- 1 + 0.5 * d$x + stats::rnorm(80)
fit <- frm(y ~ x, data = d)
nd <- data.frame(x = c(-1, 0, 1))
V <- marginaleffects::get_vcov(fit)
s1 <- marginaleffects::predictions(fit, newdata = nd)$std.error
s2 <- marginaleffects::predictions(fit, newdata = nd, vcov = 2 * V)$std.error
cat("standard errors, own:", format(s1), "\n")
cat("with vcov = 2 V:     ", format(s2), "\n")
cat("ratio:", format(s2 / s1), "(sqrt(2) = 1.414214 if the matrix is read)\n")
m <- stats::lm(y ~ x, data = d)
t1 <- marginaleffects::predictions(m, newdata = nd)$std.error
t2 <- marginaleffects::predictions(m, newdata = nd,
                                   vcov = 2 * stats::vcov(m))$std.error
cat("control, lm(): ratio", format(t2 / t1), "\n")
