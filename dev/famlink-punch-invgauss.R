# Punch item 8: `family = "inverse.gaussian"` takes brms's default link,
# 1/mu^2, where base took log. Fits the inverse gaussian designs the
# test suite already uses, with the family named by string so the
# DEFAULT link is what is exercised, in one arm:
#
#   Rscript dev/famlink-punch-invgauss.R base
#   Rscript dev/famlink-punch-invgauss.R lane
#
# Designs, each with the seed of the test it is taken from:
#   setprior   tests/testthat/test-setprior.R, seed 403, n 800, mean exp(.)
#   robust     test-numerical-robustness.R 1/mu^2 design, seed 20260916,
#              n 400, eta = 2 + 0.3 x bounded
#   band       test-numerical-robustness.R band design, seed 11, n 400
#   brmslp     test-brms-likelihood.R design, seed 20260916, n 300
arm <- commandArgs(trailingOnly = TRUE)[1]
lib <- switch(arm,
              base = "C:/Users/adf44/source/r/rellib-r3",
              lane = "C:/Users/adf44/source/r/famlink-lib",
              stop("arm must be base or lane"))
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

designs <- list(
  setprior = function() {
    set.seed(403)
    n <- 800
    x <- rnorm(n)
    data.frame(y = RTMBdist::rinvgauss(n, exp(0.5 + 0.3 * x), 2), x = x)
  },
  robust = function() {
    set.seed(20260916)
    n <- 400
    x <- rnorm(n)
    eta <- 2 + 0.3 * pmax(pmin(x, 3.5), -3.5)
    data.frame(y = stats::rgamma(n, 8, 8 * sqrt(eta)), x = x)
  },
  band = function() {
    set.seed(11)
    n <- 400
    x <- stats::rnorm(n)
    mig <- 1 / sqrt(pmax(0.5 + 0.2 * x, 0.05))
    data.frame(y = stats::rgamma(n, shape = 20, rate = 20 / mig), x = x)
  },
  brmslp = function() {
    set.seed(20260916)
    n <- 300
    x <- rnorm(n)
    eta <- 2 + 0.3 * pmax(pmin(x, 3), -3)
    data.frame(y = rgamma(n, 8, 8 * sqrt(eta)), x = x)
  }
)

cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "\n")
cat(sprintf("%-9s %-7s %-6s %-5s %-14s %-9s %s\n", "design", "link",
            "status", "conv", "logLik", "max|grad|", "warnings"))
for (nm in names(designs)) {
  d <- designs[[nm]]()
  warn <- character(0)
  fit <- withCallingHandlers(
    tryCatch(frm(y ~ x, d, family = "inverse.gaussian"),
             error = function(e) e),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (inherits(fit, "error")) {
    cat(sprintf("%-9s %-7s %-6s %s\n", nm, "", "error",
                conditionMessage(fit)))
    next
  }
  g <- tryCatch(max(abs(fit$obj$gr(fit$opt$par))), error = function(e) NA)
  cat(sprintf("%-9s %-7s %-6s %-5d %-14.6f %-9.2e %d%s\n", nm,
              family(fit)$links$mu$name, "fit", fit$opt$convergence,
              as.numeric(logLik(fit)), g, length(warn),
              if (length(warn)) paste0(" (", substr(warn[1], 1, 50), ")")
              else ""))
}
