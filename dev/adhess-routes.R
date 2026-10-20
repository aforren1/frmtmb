# Lane adhess: which covariance route the base took on the rows of
# dev/adhess-glm.R and dev/adhess-nls.R that were already exact there
# (par_units in play, or the 0.71.0 design route: a covariance factor).
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
n <- 500
show <- function(lab, f) {
  sdr <- suppressWarnings(ns$sdr_of(f))
  an <- f$cache$se_fit_an
  cat(sprintf("%-22s check tier %s | kept analysis %s | par_units %s | factor %s\n",
              lab, if (is.null(an)) "none" else an$tier,
              !is.null(f$cache$se_analysis),
              if (is.null(f$par_units)) "none" else
                paste(signif(f$par_units, 3), collapse = ","),
              !is.null(sdr$cov_fixed_root)))
}
for (k in c(6L, 10L)) for (fm in c("poisson", "bernoulli")) {
  set.seed(100 * k + match(fm, c("poisson", "bernoulli")))
  d <- if (k == 6L) data.frame(x = runif(n, 1e3, 1e3 + 1)) else
    data.frame(x = rgamma(n, shape = 2, scale = 5e3))
  X <- cbind(1, d$x)
  eta <- as.numeric(scale(d$x)) * 0.5
  d$y <- if (fm == "poisson") rpois(n, exp(0.5 + eta)) else
    rbinom(n, 1, plogis(0.2 + eta))
  f <- suppressWarnings(frm(y ~ x, family = if (fm == "poisson")
    poisson() else bernoulli(), data = d))
  show(paste0("d", k, " ", fm), f)
}
set.seed(3)
d <- data.frame(x = 20 + runif(200))
d$y <- 2 * exp(0.8 * (d$x - 20)) + rnorm(200, 0, 0.3)
f <- frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
         start = list(beta = c(2 * exp(-16), 0.8)))
show("nls x20 seed 3", f)
# whether fd_hessian()'s retake ran: optimHess() at 1e-3 has no Cholesky
# factor and some column rose by more than 1e-2
for (k in c(6L, 10L)) for (fm in c("poisson", "bernoulli")) {
  set.seed(100 * k + match(fm, c("poisson", "bernoulli")))
  d <- if (k == 6L) data.frame(x = runif(n, 1e3, 1e3 + 1)) else
    data.frame(x = rgamma(n, shape = 2, scale = 5e3))
  eta <- as.numeric(scale(d$x)) * 0.5
  d$y <- if (fm == "poisson") rpois(n, exp(0.5 + eta)) else
    rbinom(n, 1, plogis(0.2 + eta))
  f <- suppressWarnings(frm(y ~ x, family = if (fm == "poisson")
    poisson() else bernoulli(), data = d))
  at <- ns$sdr_outer_point(f$obj)
  H1 <- stats::optimHess(at, f$obj$fn, f$obj$gr)
  pdf <- !inherits(try(chol(H1), silent = TRUE), "try-error")
  rise <- 0.5 * diag(H1) * 1e-6
  cat(sprintf("d%d %-9s optimHess factors %s | rise %s\n", k, fm, pdf,
              paste(signif(rise, 3), collapse = ",")))
}
