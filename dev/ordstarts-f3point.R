# Lane ordstarts, build fix round 1: the point -404.9916 that 0.83.0's
# default starts chose on the F3 data (test-ordinal-mixture.R "a
# threshold with no standard error is said so"), for a `start =` in that
# test. Fits on rellib-r22 (0.83.0) and prints the estimates to paste.
#   Rscript dev/ordstarts-f3point.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
stopifnot(grepl("rellib-r22", find.package("frmtmb")))
set.seed(1)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- stats::rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
d$y <- 1L + rowSums(outer(lat, stats::quantile(lat, 1:3 / 4), ">"))
fit <- suppressWarnings(suppressMessages(
  frm(bf(y ~ cs(x) + z), family = mixture(sratio(), acat()), data = d)))
cat("logLik", format(as.numeric(logLik(fit)), digits = 12), "\n")
est <- fit$estimates
est$b <- NULL
dput(lapply(est, function(v) signif(v, 15)))
