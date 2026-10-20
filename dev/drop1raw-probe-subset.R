# nobs() and logLik() of update(fit) on a univariate subset() fit, with
# the data present. Arg: base | lane.
# Log: dev/drop1raw-out/probe-subset-<arm>.txt (seed 11, n = 160)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(11)
n <- 160
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1,
                sb = rep(c(TRUE, TRUE, FALSE), length.out = n))
d$z <- 1 + d$x + log(d$x2) + rnorm(n)
d$x[c(3, 6)] <- NA   # on rows the subset leaves out
fit <- q(frm(bf(z | subset(sb) ~ x + log(x2)), data = d))
u <- q(update(fit))
cat(sprintf("fit: nobs %d, rows fitted %d, logLik %.10f\n", nobs(fit),
            fit$frame$n_obs, as.numeric(logLik(fit))))
cat(sprintf("update(fit): nobs %d, rows fitted %d, logLik %.10f\n",
            nobs(u), u$frame$n_obs, as.numeric(logLik(u))))
