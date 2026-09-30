# Validation of xbeta(), zero_inflated_beta_binomial() and
# hurdle_cumulative() on the lane build. Output: dev/fams2-validate.txt.
# Seeds are set per section. Reference densities are written here from
# stats::dbeta(), stats::pbeta(), lbeta()/lchoose() and the link CDFs,
# after brms 2.23.0's Stan functions (dev/fams2-brms-src.txt), separately
# from the package code.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(RTMB)
})
say <- function(...) cat(sprintf(...), "\n", sep = "")
rel <- function(a, b) abs(a - b) / max(1, abs(b))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")

# the objective at a parameter vector, with the likelihood alone
nll_at <- function(fit, p) as.numeric(fit$obj$fn(p))
perturb <- function(fit, s) {
  set.seed(s)
  p <- fit$obj$env$last.par.best
  p + stats::rnorm(length(p), 0, 0.3)
}
at_par <- function(fit, p) {
  rnd <- fit$obj$env$random
  f2 <- fit
  f2$estimates <- fit$obj$env$parList(if (length(rnd)) p[-rnd] else p, p)
  f2
}
dpars_at <- function(fit, p) frmtmb:::eval_dpars(at_par(fit, p))[[1]]

# ---- 1. xbeta against a hand-written density --------------------------
cat("\n== 1. xbeta log-likelihood vs dbeta/pbeta ==\n")
ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi
  b <- (1 - mu) * phi
  d <- 1 + 2 * kappa
  ifelse(y <= 0, stats::pbeta(kappa / d, a, b, log.p = TRUE),
         ifelse(y >= 1,
                stats::pbeta((1 + kappa) / d, a, b, lower.tail = FALSE,
                             log.p = TRUE),
                stats::dbeta((y + kappa) / d, a, b, log = TRUE) - log(d)))
}
set.seed(101)
n <- 800
d1 <- data.frame(x = rnorm(n), g = gl(40, n / 40))
mu <- plogis(0.3 + 0.6 * d1$x)
kap <- exp(-2 + 0.5 * d1$x)
z <- rbeta(n, mu * 6, (1 - mu) * 6)
d1$y <- pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
say("  data: %d zeros, %d ones, %d interior", sum(d1$y == 0),
    sum(d1$y == 1), sum(d1$y > 0 & d1$y < 1))
f1 <- frm(bf(y ~ x, kappa ~ x), family = xbeta(), data = d1)
for (lab in c("optimum", "perturbed 1", "perturbed 2")) {
  p <- switch(lab, optimum = f1$obj$env$last.par.best,
              "perturbed 1" = perturb(f1, 1), "perturbed 2" = perturb(f1, 2))
  dp <- dpars_at(f1, p)
  r <- sum(ref_xbeta(d1$y, dp$mu, dp$phi, dp$kappa))
  o <- -nll_at(f1, p)
  say("  %-12s frmtmb %.12f  reference %.12f  rel diff %.2e", lab, o, r,
      rel(o, r))
}
for (lk in c("probit", "cloglog", "cauchit")) {
  fl <- frm(bf(y ~ x, kappa ~ x), family = xbeta(lk), data = d1)
  dp <- frmtmb:::eval_dpars(fl)[[1]]
  r <- sum(ref_xbeta(d1$y, dp$mu, dp$phi, dp$kappa))
  say("  link = %-8s at its optimum: frmtmb %.12f  reference %.12f  rel %.2e",
      lk, as.numeric(logLik(fl)), r, rel(as.numeric(logLik(fl)), r))
}
# the density alone, on the tape (with .eta_mu) and off it, over a grid
fam <- xbeta()
grid <- expand.grid(eta = c(-8, -2, 0, 2, 8), lphi = c(-1, 2, 6),
                    lkap = c(-8, -2, 1))
yv <- c(0, 1, 1e-6, 0.3, 0.97)
worst_off <- worst_on <- 0
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  mu_ <- plogis(g$eta); phi_ <- exp(g$lphi); k_ <- exp(g$lkap)
  r <- ref_xbeta(yv, mu_, phi_, k_)
  off <- fam$lpdf(yv, list(mu = mu_, phi = phi_, kappa = k_), list())
  worst_off <- max(worst_off, max(abs(off - r) / pmax(1, abs(r))))
  F <- MakeTape(function(p) {
    e <- p[1] + 0 * yv
    sum(fam$lpdf(yv, list(mu = plogis(e), .eta_mu = e, phi = exp(p[2]),
                          kappa = exp(p[3])), list()))
  }, c(0, 0, 0))
  on <- F(c(g$eta, g$lphi, g$lkap))
  worst_on <- max(worst_on, abs(on - sum(r)) / max(1, abs(sum(r))))
}
say("  grid of %d (eta, log phi, log kappa) x 5 responses: off-tape max rel %.2e, on-tape sum max rel %.2e",
    nrow(grid), worst_off, worst_on)

# ---- 2. xbeta with a random intercept: the Laplace needs third
#         derivatives, which RTMB::pbeta() does not always have ------------
cat("\n== 2. xbeta with (1 | g): the incomplete beta's third derivatives ==\n")
fam_pbeta <- frmtmb_family(
  "xbeta_pbeta", dpars = c("mu", "phi", "kappa"),
  links = list(mu = "logit", phi = "log", kappa = "log"),
  lpdf = function(y, dpars, aterms) {
    i0 <- as.numeric(y <= 0); i1 <- as.numeric(y >= 1); ib <- i0 + i1
    mp <- dpar_complement(dpars, "mu", "logit")
    a <- mp$p * dpars$phi; b <- mp$q * dpars$phi
    kap <- dpars$kappa; d <- 1 + 2 * kap
    z <- (y * (1 - ib) + 0.5 * ib + kap) / d
    out <- (1 - ib) * (RTMB::dbeta(z, a, b, log = TRUE) - log(d))
    q <- kap / d + 0 * y; a <- a + 0 * y; b <- b + 0 * y
    w0 <- which(i0 == 1); w1 <- which(i1 == 1)
    if (length(w0)) out[w0] <- log(RTMB::pbeta(q[w0], a[w0], b[w0]))
    if (length(w1)) out[w1] <- log(RTMB::pbeta(q[w1], b[w1], a[w1]))
    out
  },
  init_dpars = list(mu = function(y, a) 0.5, phi = function(y, a) 5,
                    kappa = function(y, a) 0.1))
set.seed(102)
nre <- 0
for (s in 1:10) {
  set.seed(1000 + s)
  dd <- data.frame(x = rnorm(600), g = gl(20, 30))
  u <- rnorm(20, 0, 0.4)[dd$g]
  mu <- plogis(0.3 + 0.6 * dd$x + u)
  z <- rbeta(600, mu * 6, (1 - mu) * 6)
  dd$y <- pmin(pmax(1.3 * z - 0.15, 0), 1)
  a <- tryCatch({frm(bf(y ~ x + (1 | g)), family = xbeta(), data = dd); "ok"},
                error = function(e) "error")
  b <- tryCatch({suppressWarnings(frm(bf(y ~ x + (1 | g)), family = fam_pbeta,
                                      data = dd)); "ok"},
                error = function(e) {
                  m <- conditionMessage(e)
                  if (grepl("NA/NaN gradient", m)) "NA/NaN gradient" else
                    substr(m, 1, 70)
                })
  say("  seed %d: xbeta() %s; the same density on RTMB::pbeta() %s",
      1000 + s, a, b)
}

# ---- 3. zero_inflated_beta_binomial -----------------------------------
cat("\n== 3. zero_inflated_beta_binomial vs lbeta/lchoose ==\n")
ref_bb <- function(y, n, mu, phi) {
  a <- mu * phi; b <- (1 - mu) * phi
  lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b)
}
ref_zibb <- function(y, n, mu, phi, zi) {
  base <- ref_bb(y, n, mu, phi)
  ifelse(y == 0, log(zi + (1 - zi) * exp(base)), log1p(-zi) + base)
}
set.seed(103)
n <- 800
d3 <- data.frame(x = rnorm(n), g = gl(40, n / 40),
                 tr = sample(4:20, n, TRUE))
mu <- plogis(-0.4 + 0.5 * d3$x)
yb <- rbinom(n, d3$tr, rbeta(n, mu * 4, (1 - mu) * 4))
d3$y <- ifelse(runif(n) < plogis(-1 + 0.4 * d3$x), 0L, yb)
f3 <- frm(bf(y | trials(tr) ~ x, zi ~ x),
          family = zero_inflated_beta_binomial(), data = d3)
for (lab in c("optimum", "perturbed 1", "perturbed 2")) {
  p <- switch(lab, optimum = f3$obj$env$last.par.best,
              "perturbed 1" = perturb(f3, 1), "perturbed 2" = perturb(f3, 2))
  dp <- dpars_at(f3, p)
  r <- sum(ref_zibb(d3$y, d3$tr, dp$mu, dp$phi, dp$zi))
  o <- -nll_at(f3, p)
  say("  %-12s frmtmb %.12f  reference %.12f  rel diff %.2e", lab, o, r,
      rel(o, r))
}

cat("\n== 4. zero_inflated_beta_binomial vs glmmTMB(betabinomial, ziformula) ==\n")
g3 <- glmmTMB::glmmTMB(cbind(y, tr - y) ~ x, ziformula = ~ x,
                       family = glmmTMB::betabinomial(), data = d3)
say("  fixed:  logLik frmtmb %.10f  glmmTMB %.10f  diff %.2e",
    as.numeric(logLik(f3)), as.numeric(logLik(g3)),
    as.numeric(logLik(f3)) - as.numeric(logLik(g3)))
fe <- fixef_by_dpar(f3)
gc <- glmmTMB::fixef(g3)
se_f <- sqrt(diag(vcov(f3)))
say("  fixed:  max |coef diff| mu %.2e  zi %.2e  log phi %.2e",
    max(abs(fe$mu - gc$cond)), max(abs(fe$zi - gc$zi)),
    abs(log(frmtmb:::eval_dpars(f3)[[1]]$phi[1]) - log(glmmTMB::sigma(g3))))
f3r <- frm(bf(y | trials(tr) ~ x + (1 | g), zi ~ x),
           family = zero_inflated_beta_binomial(), data = d3)
g3r <- glmmTMB::glmmTMB(cbind(y, tr - y) ~ x + (1 | g), ziformula = ~ x,
                        family = glmmTMB::betabinomial(), data = d3)
say("  (1|g):  logLik frmtmb %.10f  glmmTMB %.10f  diff %.2e",
    as.numeric(logLik(f3r)), as.numeric(logLik(g3r)),
    as.numeric(logLik(f3r)) - as.numeric(logLik(g3r)))
say("  (1|g):  max |coef diff| mu %.2e  zi %.2e",
    max(abs(fixef_by_dpar(f3r)$mu - glmmTMB::fixef(g3r)$cond)),
    max(abs(fixef_by_dpar(f3r)$zi - glmmTMB::fixef(g3r)$zi)))

# ---- 5. hurdle_cumulative ---------------------------------------------
cat("\n== 5. hurdle_cumulative vs hand-written density ==\n")
cdf <- list(logit = plogis, probit = pnorm, cauchit = pcauchy,
            cloglog = function(x) 1 - exp(-exp(x)))
ref_hc <- function(y, eta, hu, disc, tau, link) {
  F <- cdf[[link]]
  K1 <- length(tau)
  up <- ifelse(y >= K1 + 1, 1, F(disc * (tau[pmin(pmax(y, 1), K1)] - eta)))
  lo <- ifelse(y <= 1, 0, F(disc * (tau[pmax(pmin(y - 1, K1), 1)] - eta)))
  ifelse(y == 0, log(hu), log1p(-hu) + log(up - lo))
}
set.seed(105)
n <- 800
d5 <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(40, n / 40))
u <- rlogis(n) + 0.8 * d5$x
d5$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * d5$z), 0L,
               1L + (u > -1) + (u > 0.3) + (u > 1.5))
for (lk in names(cdf)) {
  f5 <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(lk), data = d5)
  for (lab in c("optimum", "perturbed 1")) {
    p <- if (lab == "optimum") f5$obj$env$last.par.best else perturb(f5, 3)
    f2 <- at_par(f5, p)
    dp <- frmtmb:::eval_dpars(f2)[[1]]
    tau <- frmtmb:::ord_tau_from_raw(f2$estimates$tau_raw, TRUE)
    r <- sum(ref_hc(d5$y, dp$mu, dp$hu, dp$disc, tau, lk))
    o <- -nll_at(f5, p)
    say("  %-8s %-12s frmtmb %.12f  reference %.12f  rel diff %.2e", lk,
        lab, o, r, rel(o, r))
  }
}
fd <- frm(bf(y ~ x, hu ~ z, disc ~ 0 + z), family = hurdle_cumulative(),
          data = d5)
dp <- frmtmb:::eval_dpars(fd)[[1]]
tau <- frmtmb:::ord_tau_from_raw(fd$estimates$tau_raw, TRUE)
r <- sum(ref_hc(d5$y, dp$mu, dp$hu, dp$disc, tau, "logit"))
say("  disc ~ 0 + z at its optimum: frmtmb %.12f  reference %.12f  rel %.2e",
    as.numeric(logLik(fd)), r, rel(as.numeric(logLik(fd)), r))

cat("\n== 6. hurdle_cumulative factorization identity (glm + MASS::polr) ==\n")
# with separate predictors for the hurdle and the ordinal part the two
# share no parameter, so the joint ML fit is the two separate fits
pos <- d5$y > 0
methods <- c(logit = "logistic", probit = "probit", cloglog = "cloglog",
             cauchit = "cauchit")
for (lk in names(methods)) {
  f6 <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(lk), data = d5)
  gl_ <- glm(I(y == 0) ~ z, family = binomial, data = d5)
  po <- MASS::polr(factor(y, ordered = TRUE) ~ x, data = d5[pos, ],
                   method = methods[[lk]])
  tot <- as.numeric(logLik(gl_)) + as.numeric(logLik(po))
  tau <- frmtmb:::ord_tau_from_raw(f6$estimates$tau_raw, TRUE)
  say("  %-8s logLik joint %.10f  glm + polr %.10f  residual %.2e  max |zeta diff| %.2e  |b_x diff| %.2e",
      lk, as.numeric(logLik(f6)), tot, as.numeric(logLik(f6)) - tot,
      max(abs(tau - po$zeta)),
      abs(fixef_by_dpar(f6)$mu[["x"]] - coef(po)[["x"]]))
}

cat("\n== 7. brms:::log_lik_hurdle_cumulative() at the frmtmb optimum ==\n")
for (lk in c("logit", "probit")) {
  f7 <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(lk), data = d5)
  dp <- frmtmb:::eval_dpars(f7)[[1]]
  tau <- frmtmb:::ord_tau_from_raw(f7$estimates$tau_raw, TRUE)
  prep <- structure(list(
    dpars = list(mu = matrix(dp$mu, 1L), hu = matrix(dp$hu, 1L),
                 disc = matrix(rep(dp$disc, length.out = n), 1L)),
    thres = list(thres = matrix(tau, 1L)),
    data = list(Y = d5$y), ndraws = 1L, nobs = n,
    family = brms::hurdle_cumulative(lk)), class = "brmsprep")
  ll <- vapply(seq_len(n), function(i) {
    as.numeric(brms:::log_lik_hurdle_cumulative(i, prep))
  }, 0)
  say("  %-7s sum frmtmb %.10f  brms %.10f  diff %.2e", lk,
      as.numeric(logLik(f7)), sum(ll), as.numeric(logLik(f7)) - sum(ll))
}

# ---- 8. fitted() against brms's posterior_epred formulas ---------------
cat("\n== 8. fitted() against brms's posterior_epred formulas ==\n")
# brms:::posterior_epred_xbeta(), term for term
epred_xbeta <- function(mu, phi, nu) {
  a <- mu * phi; b <- (1 - mu) * phi; d <- (1 + 2 * nu)
  q0 <- nu / d; q1 <- (1 + nu) / d
  t3 <- pbeta(q1, a, b)
  t1 <- d * mu * (pbeta(q1, a + 1, b) - pbeta(q0, a + 1, b))
  t2 <- nu * (t3 - pbeta(q0, a, b))
  1 + t1 - t2 - t3
}
dp <- frmtmb:::eval_dpars(f1)[[1]]
e1 <- epred_xbeta(dp$mu, dp$phi, dp$kappa)
say("  xbeta                       max rel diff %.2e",
    max(abs(frm_linpred(f1, type = "response") - e1) / abs(e1)))
dp <- frmtmb:::eval_dpars(f3)[[1]]
e3 <- dp$mu * d3$tr * (1 - dp$zi)
say("  zero_inflated_beta_binomial max rel diff %.2e",
    max(abs(frm_linpred(f3, type = "response") - e3) / abs(e3)))
f5 <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(), data = d5)
dp <- frmtmb:::eval_dpars(f5)[[1]]
tau <- frmtmb:::ord_tau_from_raw(f5$estimates$tau_raw, TRUE)
# brms:::posterior_epred_hurdle_cumulative(): cbind(hu, (1 - hu) *
# dcumulative()), with dcumulative() = inv_link_cumulative(disc (thres - eta))
Fm <- plogis(outer(-dp$mu, tau, "+"))
pc <- cbind(Fm, 1) - cbind(0, Fm)
e5 <- cbind(dp$hu, pc * (1 - dp$hu))
P5 <- frm_linpred(f5, type = "response")
say("  hurdle_cumulative           max abs diff %.2e (probabilities)",
    max(abs(unname(P5) - e5)))

# ---- 9. simulate() against the density --------------------------------
cat("\n== 9. simulate() against the fitted distribution (seed 109) ==\n")
set.seed(109)
S <- 400
sx <- simulate(f1, nsim = S)
yy <- unlist(sx, use.names = FALSE)
dp <- frmtmb:::eval_dpars(f1)[[1]]
a <- dp$mu * dp$phi; b <- (1 - dp$mu) * dp$phi; q <- dp$kappa / (1 + 2 * dp$kappa)
p0 <- pbeta(q, a, b); p1 <- pbeta(q, b, a)
m <- frmtmb:::xbeta_moments(dp$mu, dp$phi, dp$kappa)
zsc <- function(obs, expd, var) (obs - expd) / sqrt(var)
say("  xbeta  P(Y = 0): drawn %d expected %.1f z %.2f", sum(yy == 0),
    S * sum(p0), zsc(sum(yy == 0), S * sum(p0), S * sum(p0 * (1 - p0))))
say("  xbeta  P(Y = 1): drawn %d expected %.1f z %.2f", sum(yy == 1),
    S * sum(p1), zsc(sum(yy == 1), S * sum(p1), S * sum(p1 * (1 - p1))))
v <- m$m2 - m$m1^2
say("  xbeta  sum of Y: drawn %.2f expected %.2f z %.2f", sum(yy),
    S * sum(m$m1), zsc(sum(yy), S * sum(m$m1), S * sum(v)))
# the variance: mean of squared deviations against its expectation, with
# a Monte Carlo sd from the fourth moment estimated from the draws
dev2 <- (matrix(yy, ncol = S) - m$m1)^2
say("  xbeta  mean (Y - m)^2: drawn %.6f expected %.6f z %.2f",
    mean(dev2), mean(v), (mean(dev2) - mean(v)) / (sd(dev2) / sqrt(length(dev2))))
chisq_cat <- function(obs, expd) {
  keep <- expd >= 5
  x2 <- sum((obs[keep] - expd[keep])^2 / expd[keep])
  c(x2 = x2, df = sum(keep) - 1, p = stats::pchisq(x2, sum(keep) - 1,
                                                    lower.tail = FALSE))
}
sz <- simulate(f3, nsim = S)
yy <- as.matrix(sz)
dp <- frmtmb:::eval_dpars(f3)[[1]]
maxn <- max(d3$tr)
expd <- numeric(maxn + 1)
for (i in seq_len(n)) {
  pr <- exp(ref_zibb(0:d3$tr[i], d3$tr[i], dp$mu[i], dp$phi[i], dp$zi[i]))
  expd[seq_along(pr)] <- expd[seq_along(pr)] + S * pr
}
obs <- tabulate(yy + 1, maxn + 1)
cs <- chisq_cat(obs, expd)
say("  zibb   counts 0..%d pooled over %d draws: chi-square %.1f on %d df, p %.3f",
    maxn, S, cs[["x2"]], cs[["df"]], cs[["p"]])
sh <- simulate(f5, nsim = S)
yy <- as.matrix(sh)
P5 <- frm_linpred(f5, type = "response")
cs <- chisq_cat(tabulate(yy + 1, ncol(P5)), S * colSums(P5))
say("  hurdle categories 0..%d pooled over %d draws: chi-square %.2f on %d df, p %.3f",
    ncol(P5) - 1, S, cs[["x2"]], cs[["df"]], cs[["p"]])

cat("\n== 6b. cauchit: MASS::polr against ordinal::clm and frmtmb ==\n")
f6 <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative("cauchit"), data = d5)
gl_ <- glm(I(y == 0) ~ z, family = binomial, data = d5)
po <- MASS::polr(factor(y, ordered = TRUE) ~ x, data = d5[pos, ],
                 method = "cauchit")
cl <- ordinal::clm(factor(y, ordered = TRUE) ~ x, data = d5[pos, ],
                   link = "cauchit")
say("  ordinal part: polr logLik %.10f  clm logLik %.10f  frmtmb minus glm %.10f",
    as.numeric(logLik(po)), as.numeric(logLik(cl)),
    as.numeric(logLik(f6)) - as.numeric(logLik(gl_)))
say("  polr convergence code %s, clm max |gradient| %.2e",
    format(po$convergence), max(abs(cl$gradient)))
tau <- frmtmb:::ord_tau_from_raw(f6$estimates$tau_raw, TRUE)
say("  max |threshold diff| frmtmb vs clm %.2e, |b_x diff| %.2e",
    max(abs(tau - cl$alpha)), abs(fixef_by_dpar(f6)$mu[["x"]] - cl$beta[["x"]]))
ordll <- function(tau, b, link = pcauchy) {
  yy <- d5$y[pos]; eta <- b * d5$x[pos]
  up <- ifelse(yy == length(tau) + 1, 1, link(tau[pmin(yy, length(tau))] - eta))
  lo <- ifelse(yy == 1, 0, link(tau[pmax(yy - 1, 1)] - eta))
  sum(log(up - lo))
}
say("  reference ordinal log-likelihood at frmtmb's estimates %.10f, at clm's %.10f",
    ordll(tau, fixef_by_dpar(f6)$mu[["x"]]), ordll(cl$alpha, cl$beta[["x"]]))
hu_ll <- sum(dbinom(d5$y == 0, 1, frmtmb:::eval_dpars(f6)[[1]]$hu, log = TRUE))
say("  hurdle part at frmtmb's estimates %.10f, glm %.10f", hu_ll,
    as.numeric(logLik(gl_)))
