# Lane adhess, copied from dev/illse2-rev-fdcov.R (lane illse2 review,
# item 7): the filed defect "a fit without random
# effects keeps sdreport()'s finite-difference covariance".
# and on which model classes? Reference: the inverse of the AD Hessian
# (obj$he()) without random effects; a Richardson Jacobian of the AD
# gradient (numDeriv) with them.
#   ADHESS_ARM=lane|base Rscript dev/adhess-fdcov.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
set.seed(11)
n <- 300
mk <- list()
mk$exp_x20 <- function() {
  d <- data.frame(x = 20 + runif(200))
  d$y <- 2 * exp(0.8 * (d$x - 20)) + rnorm(200, 0, 0.3)
  list(f = bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), d = d,
       start = list(beta = c(2 * exp(-0.8 * 20), 0.8)))
}
mk$exp_x01 <- function() {
  d <- data.frame(x = runif(200))
  d$y <- 2 * exp(0.8 * d$x) + rnorm(200, 0, 0.3)
  list(f = bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), d = d,
       start = list(beta = c(2, 0.8)))
}
mk$logistic_nl <- function() {
  d <- data.frame(x = runif(150, 0, 20))
  d$y <- 10 / (1 + exp((10 - d$x) / 2)) + rnorm(150, 0, 0.5)
  list(f = bf(y ~ A / (1 + exp((xm - x) / sc)), A + xm + sc ~ 1,
              nl = TRUE), d = d,
       start = list(beta = c(9, 9, 1.5)))
}
mk$gauss_lm <- function() {
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  list(f = bf(y ~ x1 + x2), d = d)
}
mk$gauss_xbig <- function() {
  d <- data.frame(x = rnorm(n, 500, 100))
  d$y <- 1 + 0.002 * d$x + rnorm(n, 0, 0.1)
  list(f = bf(y ~ x), d = d)
}
mk$pois_x01 <- function() {
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(0.3 + 0.4 * d$x))
  list(f = bf(y ~ x) + poisson(), d = d)
}
mk$pois_xbig <- function() {
  d <- data.frame(x = rnorm(n, 500, 100))
  d$y <- rpois(n, exp(0.3 + 0.004 * (d$x - 500)))
  list(f = bf(y ~ x) + poisson(), d = d)
}
mk$bern_xbig <- function() {
  d <- data.frame(x = rnorm(n, 0, 100))
  d$y <- rbinom(n, 1, plogis(0.2 + 0.02 * d$x))
  list(f = bf(y ~ x) + bernoulli(), d = d)
}
mk$gamma_log <- function() {
  d <- data.frame(x = rnorm(n))
  d$y <- rgamma(n, shape = 3, rate = 3 / exp(0.5 + 0.3 * d$x))
  list(f = bf(y ~ x) + Gamma(link = "log"), d = d)
}
mk$student <- function() {
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rt(n, 5)
  list(f = bf(y ~ x) + student(), d = d)
}
mk$negbin <- function() {
  d <- data.frame(x = rnorm(n))
  d$y <- rnbinom(n, mu = exp(1 + 0.4 * d$x), size = 2)
  list(f = bf(y ~ x) + negbinomial(), d = d)
}
mk$gauss_sigma <- function() {
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n, 0, exp(0.3 * d$x))
  list(f = bf(y ~ x, sigma ~ x), d = d)
}
mk$pois_xbig_re <- function() {
  d <- data.frame(x = rnorm(n, 500, 100), g = factor(rep(1:20, length.out = n)))
  d$y <- rpois(n, exp(0.3 + 0.004 * (d$x - 500) + rnorm(20, 0, 0.3)[d$g]))
  list(f = bf(y ~ x + (1 | g)) + poisson(), d = d)
}
mk$pois_re <- function() {
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- rpois(n, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.3)[d$g]))
  list(f = bf(y ~ x + (1 | g)) + poisson(), d = d)
}
for (nm in names(mk)) {
  set.seed(match(nm, names(mk)))
  cs <- mk[[nm]]()
  r <- cap(if (is.null(cs$start)) frm(cs$f, data = cs$d) else
             frm(cs$f, data = cs$d, start = cs$start))
  f <- r$v
  if (is_err(f)) { cat(sprintf("%-14s ERROR %s\n", nm, f)); next }
  p <- f$opt$par
  re <- length(f$obj$env$random) > 0
  Href <- if (!re) f$obj$he(p) else {
    J <- numDeriv::jacobian(function(q) as.numeric(f$obj$gr(q)), p)
    (J + t(J)) / 2
  }
  Vr <- solve(Href)
  Vf <- suppressWarnings(vcov(f, full = TRUE))
  dg <- sqrt(diag(Vr))
  cor_err <- max(abs(Vf - Vr) / outer(dg, dg))
  se_err <- max(abs(sqrt(diag(Vf)) / dg - 1))
  ev <- eigen(cov2cor(Vr), only.values = TRUE)$values
  cat(sprintf("%-14s %s RE %-5s | code %d | kappa %.1e | cor err %.2e | SE err %.2e | units %s | warn %s\n",
              nm, arm, re, f$opt$convergence, max(ev) / min(ev), cor_err,
              se_err, if (is.null(f$par_units)) "none" else
                paste(signif(f$par_units, 2), collapse = ","),
              wtag(r$w)))
}
