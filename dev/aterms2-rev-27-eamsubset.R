# Reviewer, re-check: the compat table now says subset() "works" with
# every family, eam's included. Measured: univariate subset() on
# wiener and lba(3) against the fit on d[s, ], and a multivariate model
# with a wiener response and a gaussian one, each subset.
# Seeds 21, 22. Log: dev/aterms2-rev-log-27-eamsubset.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
q <- function(e) suppressWarnings(suppressMessages(e))
tr <- function(e) tryCatch(q(e), error = function(e) paste("ERROR:", conditionMessage(e)))
ll <- function(f) if (is.character(f)) f else format(as.numeric(logLik(f)), digits = 13)
set.seed(21)
dat <- frmtmb.eam:::ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.28)
dat$x <- rnorm(nrow(dat))
dat$s <- rep(c(TRUE, FALSE, TRUE), length.out = nrow(dat))
f <- tr(frm(bf(rt | vint(upper) + subset(s) ~ x, bias = 0.5),
            family = wiener(), data = dat))
g <- tr(frm(bf(rt | vint(upper) ~ x, bias = 0.5), family = wiener(),
            data = dat[dat$s, ]))
cat("wiener subset:", ll(f), " on d[s, ]:", ll(g), "\n")
if (!is.character(f)) {
  cat("  fitted rows", nrow(fitted(f)), " vs", nrow(fitted(g)),
      " simulate rows", nrow(simulate(f, seed = 1)), "\n")
}
set.seed(22)
dl <- frmtmb.eam:::lba_simulate(240, v = c(2.2, 1.4, 0.8), A = 0.5, k = 0.4,
                                ndt = 0.2)
dl$s <- rep(c(TRUE, TRUE, FALSE), length.out = nrow(dl))
fl <- tr(frm(bf(rt | vint(choice) + subset(s) ~ 1), family = lba(3), data = dl))
gl <- tr(frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = dl[dl$s, ]))
cat("lba subset:", ll(fl), " on d[s, ]:", ll(gl), "\n")
dat$y2 <- rnorm(nrow(dat))
dat$s2 <- !dat$s
fm <- tr(frm(bf(rt | vint(upper) + subset(s) ~ x, bias = 0.5) + wiener() +
               bf(y2 | subset(s2) ~ x) + gaussian(), data = dat))
cat("mv wiener + gaussian subset:", ll(fm), " sum of parts:",
    if (is.character(g)) g else format(as.numeric(logLik(g)) +
      as.numeric(logLik(q(frm(y2 ~ x, data = dat[dat$s2, ])))), digits = 13),
    "\n")
