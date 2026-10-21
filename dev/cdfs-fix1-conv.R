# Fix round 1, finding 6: the size of the change from the inclusive
# discrete censoring reading (through 0.72.0) to brms's, on one example:
# 200 poisson draws at lambda = 4 recorded as "6 or more" with
# cens(cc), cc = y >= 6, y = pmin(y, 6). The same data frame now means
# P(Y > 6) on the censored rows.
#   Rscript dev/cdfs-fix1-conv.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261020)
y <- stats::rpois(200, 4)
d <- data.frame(y = pmin(y, 6), cc = as.numeric(y >= 6))
fit <- frm(y | cens(cc) ~ 1, family = poisson(), data = d)
mu <- exp(unname(fixef_by_dpar(fit)$mu))
incl <- function(m) {
  sum(ifelse(d$cc == 1,
             stats::ppois(d$y - 1, m, lower.tail = FALSE, log.p = TRUE),
             stats::dpois(d$y, m, log = TRUE)))
}
o <- stats::optimize(function(g) -incl(exp(g)), c(0, 3), tol = 1e-12)
cat(sprintf("censored rows %d\n", sum(d$cc)))
cat(sprintf("this build: lambda %.4f logLik %.4f\n", mu,
            as.numeric(logLik(fit))))
cat(sprintf("inclusive reading: lambda %.4f logLik %.4f\n", exp(o$minimum),
            -o$objective))
cat(sprintf("logLik difference %.2f (%.3f per censored row), lambda %.2f percent\n",
            -o$objective - as.numeric(logLik(fit)),
            (-o$objective - as.numeric(logLik(fit))) / sum(d$cc),
            100 * (mu / exp(o$minimum) - 1)))
# the same "6 or more" written for this build
d2 <- data.frame(y = ifelse(y >= 6, 5, y), cc = as.numeric(y >= 6))
fit2 <- frm(y | cens(cc) ~ 1, family = poisson(), data = d2)
cat(sprintf("6 or more written as cens at 5: lambda %.4f logLik %.4f\n",
            exp(unname(fixef_by_dpar(fit2)$mu)), as.numeric(logLik(fit2))))
cat("END\n")
