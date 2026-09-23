# Shared setup for the wt-skewinit lane. LIB is set by the caller
# before sourcing, so the same scripts run against the base build
# (rellib-r3 only) and against the lane build (skewinit-lib first).
if (!exists("LIB")) LIB <- "C:/Users/adf44/source/r/skewinit-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(dirname(getNamespaceInfo("frmtmb", "path"))), "\n")

skew <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3

# The two 40-seed streams of the wt-drmtmb lane. The first run drew one
# extra rnorm(n) before the covariate; both streams are acceptance.
make_data <- function(seed, dead_draw) {
  set.seed(seed)
  n <- 200
  if (dead_draw) x <- rnorm(n)
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = y, xs = xs)
}

# The external oracle. sn::selm maximizes the same likelihood in the
# direct (DP) parameterization, so its logLik is comparable directly.
sn_ll <- function(dd) {
  m <- sn::selm(y ~ xs, family = "SN", data = dd)
  c(ll = as.numeric(m@logL),
    alpha = as.numeric(m@param$dp[["alpha"]]))
}

frm_sn <- function(dd, ...) {
  frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(),
      data = dd, ...)
}
