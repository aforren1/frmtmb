# Lane osa: where the 0.79.1 beta_binomial NaN came from, on RTMB and
# RTMBdist alone (no frmtmb). One beta-binomial regression of 150 rows at
# 10 trials (data seed 9), fitted, then oneStepPredict() three ways:
#   generic   oneStepGeneric, discrete = TRUE, range = c(0, Inf): the call
#             frmtmb 0.79.0 made for every count family
#   support   oneStepGeneric with discreteSupport = 0:10
#   cdf       method = "cdf", through RTMBdist 1.2.0's pbetabinom()
# and the same for RTMB's dbinom() on a binomial response.
#
#   Rscript dev/osa-rtmbdist.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(RTMBdist) })
cat("RTMB", format(packageVersion("RTMB")), " RTMBdist",
    format(packageVersion("RTMBdist")), "\n")
set.seed(9)
n <- 150
x <- rnorm(n)
pb <- plogis(-0.3 + 0.8 * x)
ybb <- rbinom(n, 10, rbeta(n, pb * 8, (1 - pb) * 8))
ybi <- rbinom(n, 10, pb)
fitm <- function(y, bb) {
  dat <- list(y = y, x = x)
  f <- function(p) {
    getAll(p, dat)
    y <- OBS(y)
    m <- plogis(b0 + b1 * x)
    if (bb) {
      -sum(dbetabinom(y, 10, m * exp(lphi), (1 - m) * exp(lphi), log = TRUE))
    } else {
      -sum(dbinom(y, 10, m, log = TRUE))
    }
  }
  obj <- MakeADFun(f, list(b0 = 0, b1 = 0, lphi = 2), silent = TRUE,
                   map = if (!bb) list(lphi = factor(NA)))
  nlminb(obj$par, obj$fn, obj$gr)
  obj
}
cat("BEGIN osa-rtmbdist\n")
for (bb in c(TRUE, FALSE)) {
  obj <- fitm(if (bb) ybb else ybi, bb)
  lab <- if (bb) "dbetabinom" else "dbinom"
  ways <- list(
    generic = list(method = "oneStepGeneric", discrete = TRUE,
                   range = c(0, Inf)),
    support = list(method = "oneStepGeneric", discrete = TRUE,
                   discreteSupport = 0:10),
    cdf = list(method = "cdf", discrete = TRUE))
  for (w in names(ways)) {
    r <- tryCatch(suppressWarnings(
      do.call(oneStepPredict, c(list(obj, trace = FALSE), ways[[w]]))$residual),
      error = function(e) e)
    cat(sprintf("%-11s %-8s %s\n", lab, w,
                if (inherits(r, "error")) paste("ERROR", conditionMessage(r))
                else sprintf("nan %d of %d", sum(!is.finite(r)), length(r))))
  }
}
cat("END osa-rtmbdist\n")
