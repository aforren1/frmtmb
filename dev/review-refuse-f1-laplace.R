# Short check, item 2: the Laplace that remains on read rows, against an exact
# maximum likelihood written here by quadrature over each missing value: no
# frmtmb density or chart in it. y | mi() ~ x, z ~ mi(y) + x. A reader of
# strength r: a strong reader pins the value, a weak one leaves the row nearly
# unread. Reports frm()'s shift from the exact maximum in the exact fit's SE.
# Usage: Rscript dev/review-refuse-f1-laplace.R <family>. Seed 11200 + family.
fam <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/review-refuse-lib.R")
fams <- c("Gamma", "exponential", "beta", "exgaussian", "student", "lognormal", "von_mises", "gaussian")
fi <- match(fam, fams)
ff <- switch(fam, beta = Beta(), Gamma = stats::Gamma(link = "log"), gaussian = gaussian(),
             get(fam, asNamespace("frmtmb"))())
# a density of y and a grid in the family's own natural chart
dexg <- function(y, mu, sig, bet) {
  # brms's mu is the mean: mu + exp part mean, so the gaussian part's mean is mu - bet
  m <- mu - bet
  -log(bet) + (m - y) / bet + sig^2 / (2 * bet^2) + pnorm((y - m) / sig - sig / bet, log.p = TRUE)
}
spec <- switch(fam,
  Gamma = list(np = 1, g = seq(-14, 7, length.out = 801), inv = exp, lj = identity,
    ld = function(y, eta, p) dgamma(y, exp(p), exp(p) / exp(eta), log = TRUE),
    rgen = function(eta) rgamma(length(eta), 1.3, 1.3 / exp(eta))),
  exponential = list(np = 0, g = seq(-14, 7, length.out = 801), inv = exp, lj = identity,
    ld = function(y, eta, p) dexp(y, exp(-eta), log = TRUE),
    rgen = function(eta) rexp(length(eta), exp(-eta))),
  lognormal = list(np = 1, g = seq(-14, 12, length.out = 801), inv = exp, lj = identity,
    ld = function(y, eta, p) dlnorm(y, eta, exp(p), log = TRUE),
    rgen = function(eta) rlnorm(length(eta), eta, 0.9)),
  beta = list(np = 1, g = seq(-30, 30, length.out = 1601), inv = plogis, lj = function(u) -abs(u) - 2 * log1p(exp(-abs(u))),
    ld = function(y, eta, p) { m <- plogis(eta); dbeta(y, m * exp(p), (1 - m) * exp(p), log = TRUE) },
    rgen = function(eta) rbeta(length(eta), plogis(eta) * 3, (1 - plogis(eta)) * 3)),
  exgaussian = list(np = 2, g = seq(-15, 40, length.out = 1401), inv = identity, lj = function(u) 0 * u,
    ld = function(y, eta, p) dexg(y, eta, exp(p[1]), exp(p[2])),
    rgen = function(eta) eta - 2 + rnorm(length(eta), 0, 0.5) + rexp(length(eta), 1 / 2)),
  student = list(np = 2, g = seq(-80, 80, length.out = 3201), inv = identity, lj = function(u) 0 * u,
    ld = function(y, eta, p) dt((y - eta) / exp(p[1]), 1 + exp(p[2]), log = TRUE) - p[1],
    rgen = function(eta) eta + rt(length(eta), 3)),
  gaussian = list(np = 1, g = seq(-20, 20, length.out = 1201), inv = identity, lj = function(u) 0 * u,
    ld = function(y, eta, p) dnorm(y, eta, exp(p), log = TRUE),
    rgen = function(eta) rnorm(length(eta), eta, 1)),
  von_mises = list(np = 1, g = seq(-pi, pi, length.out = 1001), inv = identity, lj = function(u) 0 * u,
    ld = function(y, eta, p) exp(p) * cos(y - 2 * atan(eta)) - log(2 * pi) - log(besselI(exp(p), 0)),
    rgen = function(eta) {
      m <- 2 * atan(eta); k <- 3
      vapply(m, function(mi) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + mi + pi) %% (2 * pi)) - pi }, 0)
    }))
exact <- function(d, start) {
  rows <- which(is.na(d$y)); ob <- setdiff(seq_len(nrow(d)), rows)
  vals <- spec$inv(spec$g); lj <- spec$lj(spec$g); dg <- diff(spec$g)[1]
  nll <- function(p) {
    by <- p[1:2]; dp <- if (spec$np) p[2 + seq_len(spec$np)]
    bz <- p[(3 + spec$np):(5 + spec$np)]; sz <- exp(p[6 + spec$np])
    ey <- by[1] + by[2] * d$x
    s <- sum(spec$ld(d$y[ob], ey[ob], dp)) +
      sum(dnorm(d$z[ob], bz[1] + bz[2] * d$y[ob] + bz[3] * d$x[ob], sz, log = TRUE))
    for (i in rows) {
      li <- spec$ld(vals, ey[i], dp) + lj + dnorm(d$z[i], bz[1] + bz[2] * vals + bz[3] * d$x[i], sz, log = TRUE)
      li[!is.finite(li)] <- -Inf
      mx <- max(li); s <- s + mx + log(sum(exp(li - mx)) * dg)
    }
    if (!is.finite(s)) return(1e10)
    -s
  }
  o <- nlminb(start, nll, control = list(iter.max = 2000, eval.max = 4000, rel.tol = 1e-11))
  H <- stats::optimHess(o$par, nll)
  list(par = o$par, se = sqrt(diag(solve(H))), ll = -o$objective, conv = o$convergence, fn = nll)
}
set.seed(11200 + fi)
for (cfg in list(c(200, 60, 0.7, 0.7), c(200, 60, 0.05, 1), c(200, 100, 0.7, 0.7), c(200, 15, 0.7, 0.7))) {
  n <- cfg[1]; m <- cfg[2]; bz <- cfg[3]; sz <- cfg[4]
  x <- rnorm(n)
  eta <- if (fam == "beta") 0.3 + 0.7 * x else if (fam == "von_mises") 0.3 + 0.5 * x else 0.3 + 0.5 * x
  y <- spec$rgen(eta)
  d <- data.frame(x = x, y = y)
  d$z <- 0.5 + bz * y + 0.3 * x + rnorm(n, 0, sz)
  d$y[sample(n, m)] <- NA
  a <- tr(frm(bf(y | mi() ~ x) + ff + bf(z ~ mi(y) + x) + gaussian(), data = d))
  if (is_err(a)) { pr("%-12s n=%d m=%d bz=%g frm: %s", fam, n, m, bz, msg(a, 140)); next }
  est <- a$estimates
  nm <- names(est$beta)
  by0 <- est$beta[grep("^y_", nm)]; bz0 <- est$beta[grep("^z_", nm)]
  dn <- names(est$betad)
  dy <- est$betad[grep("^y_", dn)]; sz0 <- est$betad[grep("^z_", dn)]
  p0 <- c(by0, dy, bz0[c(1, grep("miy", names(bz0)), grep("_x$", names(bz0)))], sz0)
  if (fam == "student") {
    # frm's chart is log(nu - 1) here; mine too
  }
  ex <- tr(exact(d, p0))
  if (is_err(ex)) { pr("%-12s exact: %s", fam, msg(ex, 140)); next }
  sh <- (p0 - ex$par) / ex$se
  labs <- c("y0", "yx", if (spec$np) paste0("dp", seq_len(spec$np)), "z0", "zy", "zx", "zs")
  # frm's own marginal at frm's maximum against my exact at the same point
  gapm <- -a$opt$objective - (-ex$fn(p0))
  pr("%-12s n=%d m=%3d reader slope %.2f: frm - exact loglik at frm's point %+.3f (per value %+.4f); shift in exact SE %s; exact conv %d",
     fam, n, m, bz, gapm, gapm / m, paste(sprintf("%s %+.2f", labs, sh), collapse = " "), ex$conv)
  pr("%-12s   frm [%s]; exact [%s]", fam, paste(round(p0, 3), collapse = ","), paste(round(ex$par, 3), collapse = ","))
}
