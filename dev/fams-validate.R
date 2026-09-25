# Validation of hurdle_negbinomial() and zero_one_inflated_beta() (lane
# fams, 2026-09-25). Every reference density below is written here from
# stats::dnbinom(), stats::pnbinom() and stats::dbeta(), separately from
# the package code, following brms 2.23.0's Stan functions
# (chunks/fun_hurdle_negbinomial.stan, fun_zero_one_inflated_beta.stan).
#
#   Rscript dev/fams-validate.R > dev/fams-validate.txt
#
# Library stack: the lane build first. Seeds are fixed per section.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-fams")
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")

rel <- function(a, b) abs(a - b) / max(1, abs(b))
say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---- references --------------------------------------------------------

# brms: y == 0 -> log(hu); y > 0 -> log(1 - hu) + NB(y) - log(1 - NB(0)).
# The normalizer is the NB upper tail P(Y > 0), which pnbinom() computes
# through the incomplete beta function without forming 1 - P(0).
ref_hnb <- function(y, mu, shape, hu) {
  ifelse(y == 0, log(hu),
         log1p(-hu) + dnbinom(y, size = shape, mu = mu, log = TRUE) -
           pnbinom(0, size = shape, mu = mu, lower.tail = FALSE,
                   log.p = TRUE))
}

# brms: y == 0 -> log(zoi) + log(1 - coi); y == 1 -> log(zoi) + log(coi);
# else log(1 - zoi) + beta(y | mu phi, (1 - mu) phi).
ref_zoib <- function(y, mu, phi, zoi, coi) {
  ifelse(y == 0, log(zoi) + log1p(-coi),
         ifelse(y == 1, log(zoi) + log(coi),
                log1p(-zoi) + dbeta(y, mu * phi, (1 - mu) * phi,
                                    log = TRUE)))
}

# ---- data ---------------------------------------------------------------

sim_hnb <- function(seed, n = 600, ngrp = 0, shape = 1.3) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  g <- factor(if (ngrp > 0) rep(seq_len(ngrp), length.out = n) else 1)
  u <- if (ngrp > 0) rnorm(ngrp, 0, 0.5)[g] else 0
  mu <- exp(0.6 + 0.4 * x + u)
  hu <- plogis(-0.4 + 0.7 * z)
  # positive part by inverse transform above the NB zero
  p0 <- dnbinom(0, size = shape, mu = mu)
  yp <- qnbinom(p0 + runif(n) * (1 - p0), size = shape, mu = mu)
  yp <- pmax(yp, 1)
  data.frame(y = ifelse(runif(n) < hu, 0, yp), x = x, z = z, g = g)
}

sim_zoib <- function(seed, n = 800) {
  set.seed(seed)
  x <- rnorm(n)
  mu <- plogis(-0.2 + 0.5 * x)
  phi <- 6
  zoi <- plogis(-1 + 0.6 * x)
  coi <- plogis(0.3 - 0.8 * x)
  edge <- runif(n) < zoi
  y <- ifelse(edge, as.numeric(runif(n) < coi),
              rbeta(n, mu * phi, (1 - mu) * phi))
  data.frame(y = y, x = x)
}

# ---- 1. hurdle_negbinomial: log-likelihood at shared points -------------

cat("\n== 1. hurdle_negbinomial log-likelihood vs dnbinom/pnbinom ==\n")
d1 <- sim_hnb(101)
f1 <- frm(bf(y ~ x, hu ~ z) + hurdle_negbinomial(), data = d1)
ref_ll_hnb <- function(par, d) {
  b <- par[names(par) == "beta"]
  bd <- par[names(par) == "betad"]
  mu <- exp(b[1] + b[2] * d$x)
  shape <- exp(bd[1])
  hu <- plogis(bd[2] + bd[3] * d$z)
  sum(ref_hnb(d$y, mu, shape, hu))
}
p_opt <- f1$obj$par
p_opt[] <- f1$opt$par
set.seed(1)
for (lab in c("optimum", "perturbed 1", "perturbed 2")) {
  p <- p_opt
  if (lab != "optimum") p[] <- p_opt + rnorm(length(p), 0, 0.3)
  a <- -f1$obj$fn(p)
  b <- ref_ll_hnb(p, d1)
  say("  %-12s frmtmb %.12f  reference %.12f  rel diff %.2e", lab, a, b,
      rel(a, b))
}
say("  logLik(fit) %.12f", as.numeric(logLik(f1)))

# the small-mu corner the normalizer is written for: eta = -25, where
# P(0) = 1 - 1.4e-11 and 1 - P(0) formed by subtraction keeps 5 digits
fam <- hurdle_negbinomial()
for (eta in c(-25, -8, 0, 4)) {
  for (tape in c(FALSE, TRUE)) {
    mu <- exp(eta)
    dp <- list(mu = mu, shape = 0.7, hu = 0.3)
    if (tape) dp <- c(dp, list(.eta_mu = eta, .eta_shape = log(0.7),
                               .eta_hu = qlogis(0.3)))
    yy <- c(0, 1, 2, 5)
    a <- as.numeric(fam$lpdf(yy, lapply(dp, rep, length.out = 4), list()))
    b <- ref_hnb(yy, mu, 0.7, 0.3)
    say("  eta %4d %-8s max rel diff %.2e", eta,
        if (tape) "on-tape" else "off-tape", max(rel(a, b)))
  }
}
# the same corner with brms's own expression, log1m((s/(mu+s))^s)
mu <- exp(-25)
brms_norm <- log1p(-(0.7 / (mu + 0.7))^0.7)
ref_norm <- pnbinom(0, size = 0.7, mu = mu, lower.tail = FALSE,
                    log.p = TRUE)
say(paste0("  eta -25: brms-form normalizer rel error %.2e (the ",
           "package's on-tape form is above)"), rel(brms_norm, ref_norm))

# a non-log mean link takes the other branch of the density
f1s <- frm(bf(y ~ x, hu ~ z) + hurdle_negbinomial(link = "sqrt"),
           data = d1)
ps <- f1s$opt$par
b <- ps[names(ps) == "beta"]
bd <- ps[names(ps) == "betad"]
a <- as.numeric(logLik(f1s))
r <- sum(ref_hnb(d1$y, (b[1] + b[2] * d1$x)^2, exp(bd[1]),
                 plogis(bd[2] + bd[3] * d1$z)))
say("  link = sqrt at its optimum: frmtmb %.12f  reference %.12f  rel %.2e",
    a, r, rel(a, r))

# ---- 2. hurdle_negbinomial vs glmmTMB truncated_nbinom2 -----------------

cat("\n== 2. hurdle_negbinomial vs glmmTMB(truncated_nbinom2, ziformula) ==\n")
suppressWarnings(suppressMessages(library(glmmTMB)))
g1 <- glmmTMB(y ~ x, ziformula = ~z, family = truncated_nbinom2,
              data = d1)
say("  fixed:  logLik frmtmb %.10f  glmmTMB %.10f  diff %.2e",
    as.numeric(logLik(f1)), as.numeric(logLik(g1)),
    as.numeric(logLik(f1)) - as.numeric(logLik(g1)))
fe <- fixef_by_dpar(f1)
say("  fixed:  max |coef diff| mu %.2e  hu %.2e  log shape %.2e",
    max(abs(fe$mu - fixef(g1)$cond)), max(abs(fe$hu - fixef(g1)$zi)),
    abs(f1$opt$par[names(f1$opt$par) == "betad"][1] -
          log(sigma(g1))))

d2 <- sim_hnb(202, n = 800, ngrp = 25)
t0 <- proc.time()[["elapsed"]]
f2 <- frm(bf(y ~ x + (1 | g), hu ~ z) + hurdle_negbinomial(), data = d2)
t_f2 <- proc.time()[["elapsed"]] - t0
g2 <- glmmTMB(y ~ x + (1 | g), ziformula = ~z, family = truncated_nbinom2,
              data = d2)
say("  (1|g):  logLik frmtmb %.10f  glmmTMB %.10f  diff %.2e",
    as.numeric(logLik(f2)), as.numeric(logLik(g2)),
    as.numeric(logLik(f2)) - as.numeric(logLik(g2)))
fe <- fixef_by_dpar(f2)
say("  (1|g):  max |coef diff| mu %.2e  hu %.2e",
    max(abs(fe$mu - fixef(g2)$cond)), max(abs(fe$hu - fixef(g2)$zi)))
say("  (1|g):  sd(g) frmtmb %.8f  glmmTMB %.8f",
    as.numeric(VarCorr(f2)[[1]]$sd[1]),
    attr(VarCorr(g2)$cond$g, "stddev"))
say("  (1|g):  frm() elapsed %.2f s (one run, a loaded 4-core box)", t_f2)

# ---- 3. zero_one_inflated_beta: log-likelihood at shared points ---------

cat("\n== 3. zero_one_inflated_beta log-likelihood vs dbeta ==\n")
d3 <- sim_zoib(303)
say("  data: %d zeros, %d ones, %d interior", sum(d3$y == 0),
    sum(d3$y == 1), sum(d3$y > 0 & d3$y < 1))
f3 <- frm(bf(y ~ x, zoi ~ x, coi ~ x, phi ~ 1) + zero_one_inflated_beta(),
          data = d3)
ref_ll_zoib <- function(par, d) {
  b <- par[names(par) == "beta"]
  bd <- par[names(par) == "betad"]
  nm <- names(f3$estimates$betad)
  bd <- setNames(bd, nm)
  mu <- plogis(b[1] + b[2] * d$x)
  phi <- exp(bd[["phi_(Intercept)"]])
  zoi <- plogis(bd[["zoi_(Intercept)"]] + bd[["zoi_x"]] * d$x)
  coi <- plogis(bd[["coi_(Intercept)"]] + bd[["coi_x"]] * d$x)
  sum(ref_zoib(d$y, mu, phi, zoi, coi))
}
p_opt <- f3$obj$par
p_opt[] <- f3$opt$par
set.seed(2)
for (lab in c("optimum", "perturbed 1", "perturbed 2")) {
  p <- p_opt
  if (lab != "optimum") p[] <- p_opt + rnorm(length(p), 0, 0.3)
  a <- -f3$obj$fn(p)
  b <- ref_ll_zoib(p, d3)
  say("  %-12s frmtmb %.12f  reference %.12f  rel diff %.2e", lab, a, b,
      rel(a, b))
}

# ---- 4. zero_one_inflated_beta: the factorization identity --------------
#
# IDENTITY. With separate predictors for mu/phi, zoi and coi, the
# log-likelihood is a sum of three terms that share no parameter:
#   sum_all  bernoulli(1{y in {0,1}} | zoi)
# + sum_edge bernoulli(y | coi)
# + sum_int  beta(y | mu, phi)
# so the joint ML fit IS the three separate ML fits, and its logLik is
# their sum. The residuals below are a numerical check of that algebra
# at optimizer precision, not an independent measurement.
cat("\n== 4. zero_one_inflated_beta factorization identity ==\n")
d3$edge <- as.numeric(d3$y == 0 | d3$y == 1)
r_zoi <- glm(edge ~ x, family = binomial, data = d3)
de <- d3[d3$edge == 1, ]
r_coi <- glm(y ~ x, family = binomial, data = de)
di <- d3[d3$edge == 0, ]
r_beta <- glmmTMB(y ~ x, family = beta_family(), data = di)
ll_sum <- as.numeric(logLik(r_zoi)) + as.numeric(logLik(r_coi)) +
  as.numeric(logLik(r_beta))
say("  logLik joint frmtmb %.10f  sum of glm + glm + glmmTMB beta %.10f",
    as.numeric(logLik(f3)), ll_sum)
say("  residual %.2e", as.numeric(logLik(f3)) - ll_sum)
fe <- fixef_by_dpar(f3)
say("  max |coef diff| mu %.2e  zoi %.2e  coi %.2e  log phi %.2e",
    max(abs(fe$mu - fixef(r_beta)$cond)),
    max(abs(fe$zoi - coef(r_zoi))), max(abs(fe$coi - coef(r_coi))),
    abs(fe$phi - log(sigma(r_beta))))

# ---- 5. epred against brms's posterior_epred_* formulas -----------------

cat("\n== 5. fitted() against brms's posterior_epred formulas ==\n")
dp <- frmtmb:::eval_dpars(f1)[[1]]
m_brms <- with(dp, mu / (1 - (shape / (mu + shape))^shape) * (1 - hu))
say("  hurdle_negbinomial max rel diff %.2e",
    max(abs(fitted(f1)[, "Estimate"] - m_brms) / m_brms))
dp <- frmtmb:::eval_dpars(f3)[[1]]
m_brms <- with(dp, zoi * coi + mu * (1 - zoi))
say("  zero_one_inflated_beta max rel diff %.2e",
    max(abs(fitted(f3)[, "Estimate"] - m_brms) / m_brms))

# ---- 6. frmtmb.sample: log_lik() at posterior draws ----------------------
#
# log_lik() evaluates the family density at each draw; the reference is
# the same hand-written density at the draw's natural-scale values.
cat("\n== 6. frmtmb.sample log_lik() against the reference, 50 draws ==\n")
suppressMessages(library(frmtmb.sample))
q <- function(expr) suppressWarnings(suppressMessages(expr))
dh <- d1[1:200, ]
fh <- frm(bf(y ~ x) + hurdle_negbinomial(), data = dh)
sh <- q(frm_sample(fh, chains = 1, iter = 100, seed = 11, refresh = 0))
mh <- as.matrix(sh)
llh <- log_lik(sh)
refh <- t(sapply(seq_len(nrow(mh)), function(i) {
  ref_hnb(dh$y, exp(mh[i, "b_Intercept"] + mh[i, "b_x"] * dh$x),
          mh[i, "shape"], mh[i, "hu"])
}))
say("  hurdle_negbinomial     max rel diff %.2e over %d x %d",
    max(abs(llh - refh) / pmax(1, abs(refh))), nrow(llh), ncol(llh))
dz <- d3[1:200, c("y", "x")]
fz <- frm(bf(y ~ x) + zero_one_inflated_beta(), data = dz)
sz <- q(frm_sample(fz, chains = 1, iter = 100, seed = 12, refresh = 0))
mz <- as.matrix(sz)
llz <- log_lik(sz)
refz <- t(sapply(seq_len(nrow(mz)), function(i) {
  ref_zoib(dz$y, plogis(mz[i, "b_Intercept"] + mz[i, "b_x"] * dz$x),
           mz[i, "phi"], mz[i, "zoi"], mz[i, "coi"])
}))
say("  zero_one_inflated_beta max rel diff %.2e over %d x %d",
    max(abs(llz - refz) / pmax(1, abs(refz))), nrow(llz), ncol(llz))

# ---- 7. brms's own R log_lik_*() at the frmtmb optimum ------------------
#
# brms 2.23.0 carries R versions of both densities for its log_lik()
# method. A minimal prep object (one draw, the dpars as 1 x n matrices)
# reaches them without compiling a model.
cat("\n== 7. brms:::log_lik_<family>() at the frmtmb optimum ==\n")
brms_ll <- function(fam, dp, y) {
  prep <- structure(list(
    dpars = lapply(dp, function(v) matrix(v, nrow = 1L)),
    data = list(Y = y), ndraws = 1L, nobs = length(y),
    family = brms::brmsfamily(fam)), class = "brmsprep")
  f <- get(paste0("log_lik_", fam), envir = asNamespace("brms"))
  vapply(seq_along(y), function(i) as.numeric(f(i, prep)), 0)
}
dp <- frmtmb:::eval_dpars(f1)[[1]]
b <- brms_ll("hurdle_negbinomial", dp[c("mu", "shape", "hu")], d1$y)
say("  hurdle_negbinomial     sum frmtmb %.10f  brms %.10f  diff %.2e",
    as.numeric(logLik(f1)), sum(b), as.numeric(logLik(f1)) - sum(b))
dp <- frmtmb:::eval_dpars(f3)[[1]]
b <- brms_ll("zero_one_inflated_beta", dp[c("mu", "phi", "zoi", "coi")],
             d3$y)
say("  zero_one_inflated_beta sum frmtmb %.10f  brms %.10f  diff %.2e",
    as.numeric(logLik(f3)), sum(b), as.numeric(logLik(f3)) - sum(b))
