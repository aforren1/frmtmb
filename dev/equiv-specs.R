# The minimal model of each sampled pair.
#
# One entry per pair, keyed "feature_a x feature_b" as sample.tsv names
# it. An entry is a function returning a spec:
#   bf, family, prior  text, evaluated once by frmtmb and once by brms,
#                      so both packages read ONE spelling; bf_frm,
#                      bf_brms, family_frm and family_brms where the two
#                      APIs differ (a custom family)
#   data, data2        the simulated data, from a fixed seed
#   frm_args, ctrl     extra arguments of frm() and of frmtmb_control()
#   brms_args          extra arguments of brms only (stanvars)
#   method, margs      the method feature(s) and their arguments (L2)
#   const              the stated L1 constant (default 0)
#   nobrms             brms has no spelling of this pair: why
#   decision           the record a frmtmb refusal rests on
#   skip               not measured, and why
# The data are small on purpose (n = 200, 10 groups) and every
# generator takes its own seed, so a pair replays alone.

eq_base <- function(n = 200, seed = 1, ngrp = 10) {
  set.seed(seed)
  d <- data.frame(
    x = round(rnorm(n), 3),
    z = round(rnorm(n), 3),
    # rows grouped and in time order within a group, as brms sorts them
    # for an autocorrelation term (an unsorted frame reorders its latent
    # vectors against the data, which the translator cannot follow)
    g = factor(sprintf("g%02d", sort(rep(seq_len(ngrp), length.out = n)))),
    f = factor(sample(c("a", "b", "c"), n, TRUE)),
    o = sample(0:3, n, TRUE))
  d$t <- stats::ave(seq_len(n), d$g, FUN = seq_along)
  u <- stats::rnorm(ngrp, 0, 0.5)
  # a random slope too, so a (1 + x | g) term is not fitted at a boundary
  v <- stats::rnorm(ngrp, 0, 0.3)
  d$eta <- 0.3 + 0.5 * d$x - 0.3 * d$z + u[as.integer(d$g)] +
    v[as.integer(d$g)] * d$x
  d
}

# a response for each family, from the common linear predictor
eq_y <- function(d, fam) {
  n <- nrow(d)
  eta <- d$eta
  switch(fam,
    gaussian = eta + rnorm(n),
    student = eta + rt(n, 5),
    lognormal = exp(0.5 * eta + rnorm(n, 0, 0.5)),
    shifted_lognormal = 0.2 + exp(0.5 * eta + rnorm(n, 0, 0.5)),
    asym_laplace = eta + rnorm(n),
    weibull = rweibull(n, 2, exp(0.5 * eta)),
    poisson = rpois(n, exp(eta)),
    negbinomial = rnbinom(n, mu = exp(eta), size = 2),
    geometric = rnbinom(n, mu = exp(eta), size = 1),
    binomial = rbinom(n, 10, plogis(eta)),
    bernoulli = rbinom(n, 1, plogis(eta)),
    beta_binomial = rbinom(n, 10, rbeta(n, 3 * plogis(eta),
                                        3 * (1 - plogis(eta)))),
    zero_inflated_negbinomial = ifelse(runif(n) < 0.2, 0L,
                                       rnbinom(n, mu = exp(eta), size = 2)),
    zero_inflated_poisson = ifelse(runif(n) < 0.2, 0L, rpois(n, exp(eta))),
    hurdle_poisson = ifelse(runif(n) < 0.25, 0L, 1L + rpois(n, exp(eta))),
    hurdle_gamma = ifelse(runif(n) < 0.2, 0, rgamma(n, 3, 3 / exp(eta))),
    hurdle_lognormal = ifelse(runif(n) < 0.2, 0,
                              exp(0.5 * eta + rnorm(n, 0, 0.5))),
    zero_inflated_beta = ifelse(runif(n) < 0.2, 0,
                                rbeta(n, 4 * plogis(eta),
                                      4 * (1 - plogis(eta)))),
    ordinal = ordered(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                          labels = 1:4)),
    von_mises = (eta + rnorm(n, 0, 0.5) + pi) %% (2 * pi) - pi,
    discrete_weibull = rgeom(n, plogis(-0.5 - 0.3 * d$x)),
    gen_extreme_value = eta + rnorm(n),
    stop("no generator for ", fam))
}

eq_dat <- function(fam = "gaussian", seed = 1, n = 200, ngrp = 10) {
  d <- eq_base(n = n, seed = seed, ngrp = ngrp)
  d$y <- eq_y(d, fam)
  d$ntr <- 10L
  d
}

# two responses; y2 shares part of y1's noise, so rescor has something
# to estimate
eq_mv <- function(seed, fam = "gaussian") {
  d <- eq_dat(fam, seed = seed)
  n <- nrow(d)
  e <- rnorm(n)
  d$y1 <- d$eta + e
  d$y2 <- 0.2 - 0.4 * d$x + 0.6 * e + rnorm(n, 0, 0.8)
  if (fam == "student") {
    d$y1 <- d$y1 + rt(n, 5) * 0.3
    d$y2 <- d$y2 + rt(n, 5) * 0.3
  }
  d
}

# a multi-membership design: two group memberships per row and a
# covariate per member for mmc()
eq_mm <- function(seed, fam = "gaussian") {
  d <- eq_dat(fam, seed = seed)
  n <- nrow(d)
  lv <- sprintf("m%02d", 1:8)
  d$g1 <- factor(sample(lv, n, TRUE), levels = lv)
  d$g2 <- factor(sample(lv, n, TRUE), levels = lv)
  d$x1 <- round(rnorm(n), 3)
  d$x2 <- round(rnorm(n), 3)
  u <- rnorm(8, 0, 0.6)
  shift <- 0.5 * (u[as.integer(d$g1)] + u[as.integer(d$g2)])
  if (fam == "gaussian") d$y <- d$y + shift
  if (fam == "ordinal") {
    d$y <- ordered(cut(d$eta + shift + rlogis(n),
                       c(-Inf, -0.5, 0.5, 1.5, Inf), labels = 1:4))
  }
  d
}

# a known covariance over the levels of g, for gr(g, cov = A)
eq_A <- function(d) {
  lv <- levels(d$g)
  k <- length(lv)
  A <- 0.5^abs(outer(seq_len(k), seq_len(k), "-"))
  dimnames(A) <- list(lv, lv)
  A
}

# a noisy predictor for me(): the truth is x, observed with sd sx
eq_me <- function(d) {
  d$sx <- 0.3
  d$xo <- round(d$x + rnorm(nrow(d), 0, 0.3), 3)
  d
}

# a two-component location mixture
eq_mix <- function(seed, fam = "gaussian") {
  d <- eq_base(seed = seed)
  n <- nrow(d)
  k <- runif(n) < 0.4
  d$y <- ifelse(k, -1.5, 1.5) + 0.5 * d$x + rnorm(n, 0, 0.6)
  if (fam == "shifted_lognormal") {
    d$y <- 0.1 + exp(ifelse(k, -0.8, 0.8) + 0.3 * d$x + rnorm(n, 0, 0.3))
  }
  d
}

# 20 sites on a 4 x 5 grid, ten rows each, for the spatial structures,
# which take coordinates built with num_factor() over one group
eq_sites <- function(d) {
  site <- rep(seq_len(20), length.out = nrow(d))
  d$pos <- frmtmb::num_factor((site - 1) %% 4, (site - 1) %/% 4)
  d$grp <- factor(1)
  d
}

eq_grid <- function(d, n = 7) {
  data.frame(x = seq(min(d$x), max(d$x), length.out = n), z = 0.2,
             f = factor("a", levels = levels(d$f)), o = 1L, t = 1L,
             xo = 0, sx = 0.3, ntr = 10L, expo = 1)
}
eq_grid2 <- function(d, n = 7) {
  g <- eq_grid(d, n)
  g$z <- -0.5
  g
}

# brms's vignette custom family, beta_binomial2, written once per API
eq_bb2_brms <- function() {
  brms::custom_family(
    "beta_binomial2", dpars = c("mu", "phi"), links = c("logit", "log"),
    lb = c(0, 0), ub = c(1, NA), type = "int", vars = "vint1[n]")
}
eq_bb2_stanvars <- function() {
  brms::stanvar(scode = paste(
    "real beta_binomial2_lpmf(int y, real mu, real phi, int T) {",
    "  return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);",
    "}",
    "int beta_binomial2_rng(real mu, real phi, int T) {",
    "  return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);",
    "}", sep = "\n"), block = "functions")
}
eq_bb2_frm <- function() {
  frmtmb::custom_family(
    "beta_binomial2", dpars = c("mu", "phi"),
    links = list(mu = "logit", phi = "log"),
    lpdf = function(y, dpars, aterms) {
      tr <- aterms$vint1
      a <- dpars$mu * dpars$phi
      b <- (1 - dpars$mu) * dpars$phi
      lchoose(tr, y) + lgamma(y + a) + lgamma(tr - y + b) -
        lgamma(tr + a + b) - lgamma(a) - lgamma(b) + lgamma(a + b)
    },
    required_aterms = "vint1", type = "discrete",
    post = list(mean_fn = function(dpars, aterms) {
      dpars$mu * aterms$vint1
    }))
}
eq_bb2 <- function(spec) {
  spec$family_frm <- "eq_bb2_frm()"
  spec$family_brms <- "eq_bb2_brms()"
  spec$brms_args$stanvars <- eq_bb2_stanvars()
  spec
}

# brms 2.23.0 keeps a flat Dirichlet(1) on a mo() simplex even under
# prior "", whose log density is lgamma(J) for J categories past the
# first (dev/brms-likelihood-tests.md; the one admitted L1 constant)
eq_mo_const <- function(J = 3, terms = 1) terms * lgamma(J)

# A frmtmb refusal is a divergence by decision only when NEWS, SPEC.md,
# dev/core-boundary.md or the brms-suite ledger gives a reason for it; a
# compat note that states the refusal is not enough. No sampled refusal
# met that standard, so no spec carries a `decision`.
eq_decision_gp <- paste(
  "exact gp(): the L1 gap is recorded, not scored; it is unexplained",
  "(matching the nugget removes about 11 of 94 nats on p094, review",
  "2026-10-28; dev/equiv-gp-gap.R cannot attribute the rest in double",
  "precision)")
# An exact gp() over 200 distinct positions has a kernel with many
# eigenvalues below frmtmb's nugget, and the L1 joint differs from brms's
# by tens to a thousand nats. The nugget explains only part of it
# (review 2026-10-28, finding 3), so L1 is recorded and not scored. The
# verdict comes from L4 where a flat-prior posterior converges in budget;
# where it did not (Rhat above 2, or past 2400 s) the pair records L1
# alone and is cannot measure.
eq_gp <- function(spec, sample = TRUE) {
  spec$l4 <- sample
  spec$l1_decision <- eq_decision_gp
  spec
}

eq_specs <- list()
S <- function(key, fn) eq_specs[[key]] <<- fn
nob <- function(what) paste(what, "has no brms spelling")

# =============================================================== untested

S("prior x nl", function() {
  d <- eq_dat(seed = 1)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)", data = d,
       prior = paste0("c(set_prior('normal(1, 2)', nlpar = 'b1'), ",
                      "set_prior('normal(0, 1)', nlpar = 'b2'))"))
})
S("prior x rescor", function() list(
  bf = "mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE)", data = eq_mv(2),
  prior = "set_prior('lkj(2)', class = 'rescor')"))
S("gaussian x |ID|", function() list(
  bf = "bf(y ~ x + (1 | p | g), sigma ~ 1 + (1 | p | g))",
  data = eq_dat(seed = 3)))
S("homcs x frm_curve", function() {
  d <- eq_dat(seed = 4)
  list(bf = "bf(y ~ x + homcs(f + 0 | g))", data = d,
       nobrms = nob("homcs()"), method = "frm_curve",
       margs = list(newdata = eq_grid(d)))
})
S("mvbf x mixture", function() {
  d <- eq_mix(5)
  d$y1 <- d$y
  d$y2 <- 0.3 * d$x + rnorm(nrow(d))
  list(bf = paste0("mvbf(bf(y1 ~ 1, mu1 ~ x, mu2 ~ x, family = ",
                   "mixture(gaussian(), gaussian())), bf(y2 ~ x, ",
                   "family = gaussian()), rescor = FALSE)"),
       family = NA, data = d)
})
S("prior x fitted", function() list(
  bf = "bf(y ~ x + z)", data = eq_dat(seed = 6),
  prior = "set_prior('normal(0, 0.5)', class = 'b')", method = "fitted"))
S("vint() x frm_curve", function() {
  d <- eq_dat("beta_binomial", seed = 7)
  eq_bb2(list(bf = "bf(y | vint(ntr) ~ x + z)", data = d,
              method = "frm_curve", margs = list(newdata = eq_grid(d))))
})
S("cens() x bounds", function() {
  d <- eq_dat(seed = 8)
  d$c <- ifelse(d$y > 1.5, "right", "none")
  d$y <- pmin(d$y, 1.5)
  # one slope, whose truth (0.5) is inside the bound: a bound that holds
  # puts the coefficient ON it, where Stan's unconstrain is -Inf
  list(bf = "bf(y | cens(c) ~ x)", data = d,
       prior = "set_prior('', class = 'b', lb = 0)")
})
S("cens() x confint_profile", function() {
  d <- eq_dat(seed = 9)
  d$c <- ifelse(d$y > 1.5, "right", "none")
  d$y <- pmin(d$y, 1.5)
  list(bf = "bf(y | cens(c) ~ x + z)", data = d, method = "confint_profile",
       margs = list(parm = "x", bpar = "b_x"))
})
S("cens() x disc", function() {
  d <- eq_dat("ordinal", seed = 10)
  d$yi <- as.integer(d$y)
  d$c <- ifelse(d$yi == 4L, "right", "none")
  list(bf = "bf(yi | cens(c) ~ x, disc ~ z)", family = "cumulative()",
       data = d)
})
S("cens() x equidistant", function() {
  d <- eq_dat("ordinal", seed = 11)
  d$yi <- as.integer(d$y)
  d$c <- ifelse(d$yi == 4L, "right", "none")
  list(bf = "bf(yi | cens(c) ~ x)",
       family = "cumulative(threshold = 'equidistant')", data = d)
})
S("cens() x simulate", function() {
  d <- eq_dat(seed = 12)
  d$c <- ifelse(d$y > 1.5, "right", "none")
  d$y <- pmin(d$y, 1.5)
  list(bf = "bf(y | cens(c) ~ x + z)", data = d, method = "simulate")
})
S("cs_pred() x frm_curve", function() {
  d <- eq_dat("ordinal", seed = 13)
  list(bf = "bf(y ~ x + cs(z))", family = "sratio()", data = d,
       method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("cs_pred() x nl", function() list(
  bf = "bf(y ~ b1 * x + cs(z), b1 ~ 1, nl = TRUE)", family = "sratio()",
  data = eq_dat("ordinal", seed = 14)))
S("cs_pred() x prior", function() list(
  bf = "bf(y ~ x + cs(z))", family = "sratio()",
  data = eq_dat("ordinal", seed = 15),
  prior = "set_prior('normal(0, 1)', class = 'b')"))
S("cs_pred() x simulate", function() list(
  bf = "bf(y ~ x + cs(z))", family = "sratio()",
  data = eq_dat("ordinal", seed = 16), method = "simulate"))
S("diag x frm_curve", function() {
  d <- eq_dat(seed = 17)
  list(bf = "bf(y ~ x + z + (1 + x || g))", data = d, method = "frm_curve",
       margs = list(newdata = eq_grid(d)))
})
S("diag x predict", function() list(
  bf = "bf(y ~ x + z + (1 + x || g))", data = eq_dat(seed = 18),
  method = "predict"))
S("equidistant x frm_curve", function() {
  d <- eq_dat("ordinal", seed = 19)
  list(bf = "bf(y ~ x + z)", family = "cumulative(threshold = 'equidistant')",
       data = d, method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("fitted x hypothesis_profile", function() list(
  bf = "bf(y ~ x + z)", data = eq_dat(seed = 20),
  method = c("fitted", "hypothesis_profile"),
  margs = list(fitted = list(), hypothesis_profile = list(h = "x + z = 0"))))
S("gen_extreme_value x mixture", function() list(
  bf = "bf(y ~ 1, mu1 ~ x, mu2 ~ x)",
  family = "mixture(gen_extreme_value(), gen_extreme_value())",
  data = eq_mix(21)))
S("gp_pred() x mvbf", function() {
  d <- eq_mv(22, "gaussian")
  d$y1 <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.4)
  eq_gp(list(bf = "mvbf(bf(y1 ~ gp(x)), bf(y2 ~ x), rescor = FALSE)",
             data = d), sample = FALSE)
})
S("gp_pred() x prior", function() {
  d <- eq_dat(seed = 23)
  d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.4)
  eq_gp(list(bf = "bf(y ~ gp(x))", data = d,
       prior = "set_prior('normal(0, 1)', class = 'sdgp')"),
        sample = FALSE)
})
S("homdiag x mvbf", function() list(
  bf = paste0("mvbf(bf(y1 ~ x + homdiag(1 + x | g)), bf(y2 ~ x), ",
              "rescor = FALSE)"), data = eq_mv(24),
  nobrms = nob("homdiag()")))
S("homtoep x predict", function() list(
  bf = "bf(y ~ x + homtoep(f + 0 | g))", data = eq_dat(seed = 25),
  nobrms = nob("homtoep()"), method = "predict"))
S("hsgp x mixture", function() {
  d <- eq_mix(26)
  list(bf = "bf(y ~ 1, mu1 ~ gp(x, k = 10), mu2 ~ x)",
       family = "mixture(gaussian(), gaussian())", data = d)
})
S("importance x equidistant", function() list(
  bf = "bf(y ~ x + (1 | g))", family = "cumulative(threshold = 'equidistant')",
  data = eq_dat("ordinal", seed = 27), frm_args = list(importance = 200L),
  l1_without = "importance"))
S("me() x arma()", function() {
  d <- eq_me(eq_dat(seed = 28))
  list(bf = "bf(y ~ me(xo, sx) + z + arma(time = t, gr = g))", data = d)
})
S("mi() x frm_curve", function() {
  d <- eq_dat(seed = 29)
  d$y[sample(nrow(d), 20)] <- NA
  list(bf = "bf(y | mi() ~ x + z)", data = d, method = "frm_curve",
       margs = list(newdata = eq_grid(d)))
})
S("mo() x rescor", function() {
  d <- eq_mv(30)
  d$y1 <- d$y1 + c(0, 0.5, 0.8, 1)[d$o + 1]
  list(bf = "mvbf(bf(y1 ~ mo(o) + x), bf(y2 ~ x), rescor = TRUE)",
       data = d, const = eq_mo_const())
})
S("nl x hypothesis_profile", function() {
  d <- eq_dat(seed = 31)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)", data = d,
       method = "hypothesis_profile",
       margs = list(h = "b2_Intercept = 0.3"))
})
S("nl x simulate", function() {
  d <- eq_dat(seed = 32)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)", data = d,
       method = "simulate")
})
S("ou x mixture", function() list(
  bf = "bf(y ~ 1 + ou(f + 0 | g), mu1 ~ x, mu2 ~ x)",
  family = "mixture(gaussian(), gaussian())", data = eq_mix(33),
  nobrms = nob("ou()")))
S("predict x residuals", function() list(
  bf = "bf(y ~ x + z)", data = eq_dat(seed = 34),
  method = c("predict", "residuals")))
S("prior x emmeans", function() list(
  bf = "bf(y ~ f + x)", data = eq_dat(seed = 35),
  prior = "set_prior('normal(0, 0.5)', class = 'b')", method = "emmeans",
  margs = list(specs = ~f)))
S("prior x frm_curve", function() {
  d <- eq_dat(seed = 36)
  list(bf = "bf(y ~ x + z)", data = d,
       prior = "set_prior('normal(0, 0.3)', class = 'b')",
       method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("prior x hypothesis_profile", function() list(
  bf = "bf(y ~ x + z)", data = eq_dat(seed = 37),
  prior = "set_prior('normal(0, 0.3)', class = 'b')",
  method = "hypothesis_profile", margs = list(h = "x = 0")))
S("prior x mvbf", function() list(
  bf = "mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = FALSE)", data = eq_mv(38),
  prior = "set_prior('normal(0, 0.3)', class = 'b', resp = 'y1')"))
S("profile x frm_curve", function() {
  d <- eq_dat(seed = 39)
  list(bf = "bf(y ~ x + z + (1 | g))", data = d,
       ctrl = list(profile = TRUE), method = "frm_curve",
       margs = list(newdata = eq_grid(d)))
})
S("rate() x residuals", function() {
  d <- eq_dat(seed = 40)
  d$expo <- round(runif(nrow(d), 0.5, 3), 2)
  d$y <- rpois(nrow(d), d$expo * exp(d$eta))
  list(bf = "bf(y | rate(expo) ~ x + z)", family = "poisson()", data = d,
       method = "residuals")
})
S("rescor x mm()", function() {
  d <- eq_mm(41)
  d$y1 <- d$y
  d$y2 <- 0.5 * d$y + rnorm(nrow(d))
  list(bf = "mvbf(bf(y1 ~ x + (1 | mm(g1, g2))), bf(y2 ~ x), rescor = TRUE)",
       data = d)
})
S("s() x nl", function() {
  d <- eq_dat(seed = 42)
  d$y <- (2 + sin(d$z)) * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 ~ s(z), b2 ~ 1, nl = TRUE)",
       data = d)
})
S("s() x prior", function() {
  d <- eq_dat(seed = 43)
  d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ s(x))", data = d,
       prior = "set_prior('student_t(3, 0, 1)', class = 'sds')")
})
S("s() x REML", function() {
  d <- eq_dat(seed = 44)
  d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ s(x))", data = d, frm_args = list(REML = TRUE))
})
S("shifted_lognormal x mixture", function() list(
  bf = "bf(y ~ 1, mu1 ~ x, mu2 ~ x)",
  family = "mixture(shifted_lognormal(), shifted_lognormal())",
  data = eq_mix(45, "shifted_lognormal")))
S("simulate x emmeans", function() list(
  bf = "bf(y ~ f + x)", data = eq_dat(seed = 46),
  method = c("simulate", "emmeans"),
  margs = list(simulate = list(), emmeans = list(specs = ~f))))
S("t2() x nl", function() {
  d <- eq_dat(seed = 47)
  d$y <- (2 + sin(d$z) * cos(d$x)) * exp(0.3 * d$x) +
    rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 ~ t2(z, x), b2 ~ 1, nl = TRUE)",
       data = d)
})
S("t2() x rescor", function() {
  d <- eq_mv(48)
  d$y1 <- d$y1 + sin(d$z) * cos(d$x)
  list(bf = "mvbf(bf(y1 ~ t2(x, z)), bf(y2 ~ x), rescor = TRUE)", data = d)
})
S("thres() x nl", function() list(
  bf = "bf(y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE)",
  family = "cumulative()", data = eq_dat("ordinal", seed = 49)))
S("trials() x prior", function() list(
  bf = "bf(y | trials(ntr) ~ x + z)", family = "binomial()",
  data = eq_dat("binomial", seed = 50),
  prior = "set_prior('normal(0, 0.5)', class = 'b')"))
S("trials() x ps()", function() list(
  bf = "bf(y | trials(ntr) ~ ps(x, k = 8))", family = "binomial()",
  data = eq_dat("binomial", seed = 51), nobrms = nob("ps()")))
S("trials() x simulate", function() list(
  bf = "bf(y | trials(ntr) ~ x + z)", family = "binomial()",
  data = eq_dat("binomial", seed = 52), method = "simulate"))
S("us x emmeans", function() list(
  bf = "bf(y ~ f + x + (1 + x | g))", data = eq_dat(seed = 53),
  method = "emmeans", margs = list(specs = ~f)))
S("us x frm_curve", function() {
  d <- eq_dat(seed = 54)
  list(bf = "bf(y ~ x + z + (1 + x | g))", data = d, method = "frm_curve",
       margs = list(newdata = eq_grid(d)))
})
S("us x mixture", function() list(
  bf = "bf(y ~ 1, mu1 ~ x + (1 | g), mu2 ~ x)",
  family = "mixture(gaussian(), gaussian())", data = eq_mix(55)))
S("us x predict", function() list(
  bf = "bf(y ~ x + z + (1 + x | g))", data = eq_dat(seed = 56),
  method = "predict"))
S("us x simulate", function() list(
  bf = "bf(y ~ x + z + (1 + x | g))", family = "poisson()",
  data = eq_dat("poisson", seed = 57), method = "simulate"))
S("vint() x prior", function() eq_bb2(list(
  bf = "bf(y | vint(ntr) ~ x + z)", data = eq_dat("beta_binomial", seed = 58),
  prior = "set_prior('normal(0, 0.5)', class = 'b')")))
S("vint() x sum_to_zero", function() list(
  bf = "bf(y | vint(ntr) ~ x)",
  family = "cumulative(threshold = 'sum_to_zero')",
  data = eq_dat("ordinal", seed = 59)))
S("von_mises x |ID|", function() list(
  bf = "bf(y ~ x + (1 | p | g), kappa ~ 1 + (1 | p | g))",
  family = "von_mises()", data = eq_dat("von_mises", seed = 60)))

# ============================================================ conditional

S("gaussian x hypothesis_profile", function() list(
  bf = "bf(y ~ x + z)", data = eq_dat(seed = 61),
  method = "hypothesis_profile", margs = list(h = "x - z = 0")))
S("gaussian x smooth", function() {
  d <- eq_dat(seed = 62)
  d$y <- sin(2 * d$x) + 0.3 * d$z + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ s(x) + z)", data = d)
})
S("us x smooth", function() {
  d <- eq_dat(seed = 63)
  d$y <- d$y + sin(2 * d$z)
  list(bf = "bf(y ~ s(z) + x + (1 + x | g))", data = d)
})
S("cumulative x frm_curve", function() {
  d <- eq_dat("ordinal", seed = 64)
  list(bf = "bf(y ~ x + z)", family = "cumulative()", data = d,
       method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("gaussian x mo()", function() {
  d <- eq_dat(seed = 65)
  d$y <- d$y + c(0, 0.5, 0.8, 1)[d$o + 1]
  list(bf = "bf(y ~ mo(o) + x)", data = d, const = eq_mo_const())
})
S("poisson x frm_curve", function() {
  d <- eq_dat("poisson", seed = 66)
  list(bf = "bf(y ~ x + z)", family = "poisson()", data = d,
       method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("smooth x simulate", function() {
  d <- eq_dat(seed = 67)
  d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ s(x))", data = d, method = "simulate")
})
S("student x rescor", function() list(
  bf = "mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE)", family = "student()",
  data = eq_mv(68, "student")))
S("us x importance", function() list(
  bf = "bf(y ~ x + (1 + x | g))", family = "poisson()",
  data = eq_dat("poisson", seed = 69), frm_args = list(importance = 200L),
  l1_without = "importance"))
S("|ID| x mmc()", function() list(
  bf = paste0("bf(y ~ x + (1 + mmc(x1, x2) | p | mm(g1, g2)), ",
              "sigma ~ 1 + (1 | p | mm(g1, g2)))"), data = eq_mm(70)))
S("ar1 x prior", function() list(
  bf = "bf(y ~ x + ar1(f + 0 | g))", data = eq_dat(seed = 71),
  prior = "set_prior('normal(0, 1)', class = 'b')", nobrms = nob("ar1()")))
S("asym_laplace x gr_cov", function() {
  d <- eq_dat("asym_laplace", seed = 72)
  list(bf = "bf(y ~ x + (1 | gr(g, cov = A)))", family = "asym_laplace()",
       data = d, data2 = list(A = eq_A(d)))
})
S("bernoulli x fitted", function() list(
  bf = "bf(y ~ x + z)", family = "bernoulli()",
  data = eq_dat("bernoulli", seed = 73), method = "fitted"))
S("binomial x residuals", function() list(
  bf = "bf(y | trials(ntr) ~ x + z)", family = "binomial()",
  data = eq_dat("binomial", seed = 74), method = "residuals"))
S("cumulative x hypothesis_profile", function() list(
  bf = "bf(y ~ x + z)", family = "cumulative()",
  data = eq_dat("ordinal", seed = 75), method = "hypothesis_profile",
  margs = list(h = "x = 0")))
S("cumulative x smooth", function() {
  d <- eq_dat("ordinal", seed = 76)
  list(bf = "bf(y ~ s(x) + z)", family = "cumulative()", data = d)
})
S("discrete_weibull x emmeans", function() list(
  bf = "bf(y ~ f + x)", family = "discrete_weibull()",
  data = eq_dat("discrete_weibull", seed = 77), method = "emmeans",
  margs = list(specs = ~f)))
S("gau x mo()", function() list(
  bf = "bf(y ~ mo(o) + gau(pos + 0 | grp))",
  data = eq_sites(eq_dat(seed = 78)), nobrms = nob("gau()")))
S("gaussian x fitted", function() list(
  bf = "bf(y ~ x + z, sigma ~ x)", data = eq_dat(seed = 79),
  method = "fitted"))
S("gaussian x frm_curve", function() {
  d <- eq_dat(seed = 80)
  list(bf = "bf(y ~ x + z)", data = d, method = "frm_curve",
       margs = list(newdata = eq_grid(d)))
})
S("gaussian x gr_cov", function() {
  d <- eq_dat(seed = 81)
  list(bf = "bf(y ~ x + (1 | gr(g, cov = A)))", data = d,
       data2 = list(A = eq_A(d)))
})
S("gaussian x hsgp", function() {
  d <- eq_dat(seed = 82)
  d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.4)
  list(bf = "bf(y ~ gp(x, k = 10))", data = d)
})
S("gaussian x mat", function() list(
  bf = "bf(y ~ x + mat(pos + 0 | grp))",
  data = eq_sites(eq_dat(seed = 83)), nobrms = nob("mat()")))
S("gaussian x mmc()", function() list(
  bf = "bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2)))", data = eq_mm(84)))
S("geometric x confint_profile", function() list(
  bf = "bf(y ~ x + z)", family = "geometric()",
  data = eq_dat("geometric", seed = 85), method = "confint_profile",
  margs = list(parm = "x", bpar = "b_x")))
S("gr_cov x fitted", function() {
  d <- eq_dat(seed = 86)
  list(bf = "bf(y ~ x + (1 | gr(g, cov = A)))", data = d,
       data2 = list(A = eq_A(d)), method = "fitted")
})
S("gr_cov x mmc()", function() {
  d <- eq_mm(87)
  lv <- levels(d$g1)
  A <- 0.5^abs(outer(seq_along(lv), seq_along(lv), "-"))
  dimnames(A) <- list(lv, lv)
  list(bf = "bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2, cov = A)))",
       data = d, data2 = list(A = A))
})
S("hetar1 x me()", function() list(
  bf = "bf(y ~ me(xo, sx) + hetar1(f + 0 | g))",
  data = eq_me(eq_dat(seed = 88)), nobrms = nob("hetar1()")))
S("hetar1 x prior", function() list(
  bf = "bf(y ~ x + hetar1(f + 0 | g))", data = eq_dat(seed = 89),
  prior = "set_prior('normal(0, 1)', class = 'b')",
  nobrms = nob("hetar1()")))
S("homcs x smooth", function() list(
  bf = "bf(y ~ s(x) + homcs(f + 0 | g))", data = eq_dat(seed = 90),
  nobrms = nob("homcs()")))
S("homdiag x equalto", function() {
  V <- diag(c(0.3, 0.3, 0.3))
  list(bf = "bf(y ~ x + homdiag(1 | g) + equalto(f + 0 | g, V))",
       data = eq_dat(seed = 91), data2 = list(V = V),
       nobrms = nob("homdiag() and equalto()"))
})
S("hsgp x fitted", function() {
  d <- eq_dat(seed = 92)
  d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.4)
  list(bf = "bf(y ~ gp(x, k = 10))", data = d, method = "fitted")
})
S("hurdle_cumulative x residuals", function() {
  d <- eq_dat("ordinal", seed = 93)
  d$y <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)),
                 levels = 0:4)
  list(bf = "bf(y ~ x + z)", family = "hurdle_cumulative()", data = d,
       method = "residuals")
})
S("lognormal x gp", function() {
  d <- eq_dat(seed = 94)
  d$y <- exp(0.5 * sin(2 * d$x) + rnorm(nrow(d), 0, 0.3))
  eq_gp(list(bf = "bf(y ~ gp(x))", family = "lognormal()", data = d),
        sample = FALSE)
})
S("me() x prior", function() list(
  bf = "bf(y ~ me(xo, sx) + z)", data = eq_me(eq_dat(seed = 95)),
  prior = "set_prior('normal(0, 1)', class = 'b', coef = 'mexosx')"))
S("negbinomial x fitted", function() list(
  bf = "bf(y ~ x + z)", family = "negbinomial()",
  data = eq_dat("negbinomial", seed = 96), method = "fitted"))
S("negbinomial x frm_curve_contrast", function() {
  d <- eq_dat("negbinomial", seed = 97)
  list(bf = "bf(y ~ x + z)", family = "negbinomial()", data = d,
       method = "frm_curve_contrast",
       margs = list(newdata = eq_grid(d), contrast = eq_grid2(d)))
})
S("negbinomial x smooth", function() {
  d <- eq_dat(seed = 98)
  d$y <- rnbinom(nrow(d), mu = exp(0.5 + sin(2 * d$x)), size = 3)
  list(bf = "bf(y ~ s(x))", family = "negbinomial()", data = d)
})
S("poisson x emmeans", function() list(
  bf = "bf(y ~ f + x)", family = "poisson()",
  data = eq_dat("poisson", seed = 99), method = "emmeans",
  margs = list(specs = ~f)))
S("poisson x mo()", function() {
  d <- eq_dat(seed = 100)
  d$y <- rpois(nrow(d), exp(d$eta + c(0, 0.3, 0.5, 0.6)[d$o + 1]))
  list(bf = "bf(y ~ mo(o) + x)", family = "poisson()", data = d,
       const = eq_mo_const())
})
S("poisson x residuals", function() list(
  bf = "bf(y ~ x + z)", family = "poisson()",
  data = eq_dat("poisson", seed = 101), method = "residuals"))
S("predict x mmc()", function() list(
  bf = "bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2)))", data = eq_mm(102),
  method = "predict"))
S("prior x mmc()", function() list(
  bf = "bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2)))", data = eq_mm(103),
  prior = "set_prior('normal(0, 0.5)', class = 'sd')"))
S("shifted_lognormal x frm_curve", function() {
  d <- eq_dat("shifted_lognormal", seed = 104)
  list(bf = "bf(y ~ x + z)", family = "shifted_lognormal()", data = d,
       method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("smooth x hypothesis_profile", function() {
  d <- eq_dat(seed = 105)
  d$y <- sin(2 * d$x) + 0.3 * d$z + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ s(x) + z)", data = d, method = "hypothesis_profile",
       margs = list(h = "z = 0.3"))
})
S("smooth x mmc()", function() {
  d <- eq_mm(106)
  d$y <- d$y + sin(2 * d$z)
  list(bf = "bf(y ~ s(z) + x + (1 + mmc(x1, x2) | mm(g1, g2)))", data = d)
})
S("sratio x frm_curve", function() {
  d <- eq_dat("ordinal", seed = 107)
  list(bf = "bf(y ~ x + z)", family = "sratio()", data = d,
       method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("student x frm_curve", function() {
  d <- eq_dat("student", seed = 108)
  list(bf = "bf(y ~ x + z)", family = "student()", data = d,
       method = "frm_curve", margs = list(newdata = eq_grid(d)))
})
S("student x gp", function() {
  d <- eq_dat(seed = 109)
  d$y <- sin(2 * d$x) + rt(nrow(d), 5) * 0.4
  eq_gp(list(bf = "bf(y ~ gp(x))", family = "student()", data = d),
        sample = FALSE)
})
S("student x smooth", function() {
  d <- eq_dat(seed = 110)
  d$y <- sin(2 * d$x) + rt(nrow(d), 5) * 0.4
  list(bf = "bf(y ~ s(x))", family = "student()", data = d)
})
S("sum_to_zero x mmc()", function() list(
  bf = "bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2)))",
  family = "cumulative(threshold = 'sum_to_zero')",
  data = eq_mm(111, "ordinal")))
S("trunc() x ma()", function() {
  d <- eq_dat(seed = 112)
  d <- d[d$y > -1.5, ]
  list(bf = "bf(y | trunc(lb = -1.5) ~ x + z + ma(time = t, gr = g))",
       data = d)
})
S("us x gr_cov", function() {
  d <- eq_dat(seed = 113)
  list(bf = "bf(y ~ x + (1 + x | gr(g, cov = A)))", data = d,
       data2 = list(A = eq_A(d)))
})
S("us x hsgp", function() {
  d <- eq_dat(seed = 114)
  d$y <- d$y + sin(2 * d$z)
  list(bf = "bf(y ~ gp(z, k = 10) + x + (1 + x | g))", data = d)
})
S("us x me()", function() list(
  bf = "bf(y ~ me(xo, sx) + z + (1 + z | g))",
  data = eq_me(eq_dat(seed = 115))))
S("us x mmc()", function() list(
  bf = "bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2)))", family = "poisson()",
  data = {
    d <- eq_mm(116)
    d$y <- rpois(nrow(d), exp(0.3 + 0.4 * d$x + 0.2 * (d$x1 + d$x2) / 2))
    d
  }))
S("zero_inflated_negbinomial x hsgp", function() {
  d <- eq_dat(seed = 117)
  d$y <- ifelse(runif(nrow(d)) < 0.2, 0L,
                rnbinom(nrow(d), mu = exp(0.5 + sin(2 * d$x)), size = 3))
  list(bf = "bf(y ~ gp(x, k = 10))", family = "zero_inflated_negbinomial()",
       data = d)
})

# ================================================================== works

S("gaussian x cens()", function() {
  d <- eq_dat(seed = 118)
  d$c <- ifelse(d$y > 1.5, 1L, ifelse(d$y < -1, -1L, 0L))
  d$y <- pmax(pmin(d$y, 1.5), -1)
  list(bf = "bf(y | cens(c) ~ x + z)", data = d)
})
S("gaussian x nl", function() {
  d <- eq_dat(seed = 119)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)", data = d)
})
S("gaussian x us", function() list(
  bf = "bf(y ~ x + (1 + x | g))", data = eq_dat(seed = 120)))
S("poisson x double_bar", function() list(
  bf = "bf(y ~ x + (1 + x || g))", family = "poisson()",
  data = eq_dat("poisson", seed = 121)))
S("us x cens()", function() {
  d <- eq_dat(seed = 122)
  d$c <- ifelse(d$y > 1.5, "right", "none")
  d$y <- pmin(d$y, 1.5)
  list(bf = "bf(y | cens(c) ~ x + (1 + x | g))", data = d)
})
S("beta_binomial x verbose", function() list(
  bf = "bf(y | trials(ntr) ~ x + z)", family = "beta_binomial()",
  data = eq_dat("beta_binomial", seed = 123),
  frm_args = list(verbose = TRUE)))
S("binomial x importance", function() list(
  bf = "bf(y | trials(ntr) ~ x + (1 | g))", family = "binomial()",
  data = eq_dat("binomial", seed = 124), frm_args = list(importance = 200L),
  l1_without = "importance"))
S("binomial x me()", function() {
  d <- eq_me(eq_dat("binomial", seed = 125))
  list(bf = "bf(y | trials(ntr) ~ me(xo, sx) + z)", family = "binomial()",
       data = d)
})
S("cens() x index()", function() {
  d <- eq_dat(seed = 126)
  d$c <- ifelse(d$y > 1.5, "right", "none")
  d$y <- pmin(d$y, 1.5)
  d$id <- seq_len(nrow(d))
  list(bf = "bf(y | cens(c) + index(id) ~ x + z)", data = d)
})
S("cens() x nl", function() {
  d <- eq_dat(seed = 127)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  d$c <- ifelse(d$y > 3.5, "right", "none")
  d$y <- pmin(d$y, 3.5)
  list(bf = "bf(y | cens(c) ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)",
       data = d)
})
S("cosy() x mvbf", function() list(
  bf = "mvbf(bf(y1 ~ x + cosy(gr = g)), bf(y2 ~ x), rescor = FALSE)",
  data = eq_mv(128)))
S("cratio x mvbf", function() {
  d <- eq_mv(129)
  d$o1 <- eq_y(d, "ordinal")
  list(bf = paste0("mvbf(bf(o1 ~ x, family = cratio()), ",
                   "bf(y2 ~ x, family = gaussian()), rescor = FALSE)"),
       family = NA, data = d)
})
S("cs x se()", function() {
  d <- eq_dat(seed = 130)
  d$s <- 0.5
  list(bf = "bf(y | se(s, sigma = TRUE) ~ x + cs(f + 0 | g))", data = d,
       nobrms = nob("the cs() covariance structure"))
})
S("cs_pred() x mixture", function() list(
  bf = "bf(y ~ 1, mu1 ~ cs(x), mu2 ~ cs(x))",
  family = "mixture(sratio(), sratio())",
  data = eq_dat("ordinal", seed = 131)))
S("gaussian x cosy()", function() list(
  bf = "bf(y ~ x + z + cosy(gr = g))", data = eq_dat(seed = 132)))
S("gaussian x diag", function() list(
  bf = "bf(y ~ x + (1 + x || g))", data = eq_dat(seed = 133)))
S("gaussian x prior", function() list(
  bf = "bf(y ~ x + z)", data = eq_dat(seed = 134),
  prior = paste0("c(set_prior('normal(0, 0.5)', class = 'b'), ",
                 "set_prior('student_t(3, 0, 2)', class = 'sigma'))")))
S("gaussian x simulate", function() list(
  bf = "bf(y ~ x + z)", data = eq_dat(seed = 135), method = "simulate"))
S("gp x nl", function() {
  d <- eq_dat(seed = 136)
  # z on 0.1 steps (about 40 positions): over 200 distinct positions the
  # pair did not finish inside the 2400 s budget of the first run
  d$z <- round(d$z, 1)
  d$y <- (2 + sin(d$z)) * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  eq_gp(list(bf = "bf(y ~ b1 * exp(b2 * x), b1 ~ gp(z), b2 ~ 1, nl = TRUE)",
       data = d))
})
S("lognormal x simulate", function() list(
  bf = "bf(y ~ x + z)", family = "lognormal()",
  data = eq_dat("lognormal", seed = 137), method = "simulate"))
S("mo() x mm()", function() {
  d <- eq_mm(138)
  d$y <- d$y + c(0, 0.5, 0.8, 1)[d$o + 1]
  list(bf = "bf(y ~ mo(o) + x + (1 | mm(g1, g2)))", data = d,
       const = eq_mo_const())
})
S("multinomial x us", function() {
  d <- eq_base(seed = 139)
  n <- nrow(d)
  p <- cbind(1, exp(0.3 + 0.5 * d$x), exp(-0.2 + d$eta))
  p <- p / rowSums(p)
  d$ntr <- 10L
  d$y <- t(vapply(seq_len(n), function(i) {
    stats::rmultinom(1, 10, p[i, ])[, 1]
  }, numeric(3)))
  colnames(d$y) <- c("a", "b", "c")
  list(bf = "bf(y | trials(ntr) ~ x + (1 | g))", family = "multinomial()",
       data = d)
})
S("poisson x prior", function() list(
  bf = "bf(y ~ x + z)", family = "poisson()",
  data = eq_dat("poisson", seed = 140),
  prior = "set_prior('normal(0, 0.3)', class = 'b')"))
S("poisson x simulate", function() list(
  bf = "bf(y ~ x + z)", family = "poisson()",
  data = eq_dat("poisson", seed = 141), method = "simulate"))
S("rr x frm_curve", function() {
  d <- eq_dat(seed = 142)
  list(bf = "bf(y ~ x + z + rr(f + 0 | g, d = 1))", data = d,
       nobrms = nob("rr()"), method = "frm_curve",
       margs = list(newdata = eq_grid(d)))
})
S("smooth x nl", function() {
  d <- eq_dat(seed = 143)
  d$y <- (2 + sin(d$z)) * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 ~ s(z), b2 ~ 1, nl = TRUE)",
       data = d)
})
S("sratio x prior", function() list(
  bf = "bf(y ~ x + z)", family = "sratio()",
  data = eq_dat("ordinal", seed = 144),
  prior = paste0("c(set_prior('normal(0, 0.5)', class = 'b'), ",
                 "set_prior('normal(0, 2)', class = 'Intercept'))")))
S("thres() x prior", function() list(
  bf = "bf(y | thres(3) ~ x + z)", family = "cumulative()",
  data = eq_dat("ordinal", seed = 145),
  prior = "set_prior('normal(0, 2)', class = 'Intercept')"))
S("trials() x cs_pred()", function() list(
  bf = "bf(y | trials(ntr) ~ x + cs(z))", family = "binomial()",
  data = eq_dat("binomial", seed = 146)))
S("trials() x mi_pred()", function() {
  d <- eq_dat("binomial", seed = 147)
  d$xm <- d$x
  d$xm[sample(nrow(d), 20)] <- NA
  list(bf = paste0("mvbf(bf(y | trials(ntr) ~ mi(xm) + z, ",
                   "family = binomial()), bf(xm | mi() ~ z, ",
                   "family = gaussian()), rescor = FALSE)"),
       family = NA, data = d)
})
S("trunc() x s()", function() {
  d <- eq_dat(seed = 148)
  d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.5)
  d <- d[d$y > -1, ]
  list(bf = "bf(y | trunc(lb = -1) ~ s(x))", data = d)
})
S("us x mo()", function() {
  d <- eq_dat(seed = 149)
  d$y <- d$y + c(0, 0.5, 0.8, 1)[d$o + 1]
  list(bf = "bf(y ~ mo(o) + x + (1 + x | g))", data = d,
       const = eq_mo_const())
})
S("us x trials()", function() list(
  bf = "bf(y | trials(ntr) ~ x + (1 + x | g))", family = "binomial()",
  data = eq_dat("binomial", seed = 150)))
S("verbose x rescor", function() list(
  bf = "mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE)", data = eq_mv(151),
  frm_args = list(verbose = TRUE)))
S("weibull x cens()", function() {
  d <- eq_dat("weibull", seed = 152)
  d$c <- ifelse(d$y > 2, "right", "none")
  d$y <- pmin(d$y, 2)
  list(bf = "bf(y | cens(c) ~ x + z)", family = "weibull()", data = d)
})

# ================================================================ refused

S("acat x rescor", function() {
  d <- eq_mv(153)
  d$o1 <- eq_y(d, "ordinal")
  list(bf = paste0("mvbf(bf(o1 ~ x, family = acat()), ",
                   "bf(y2 ~ x, family = gaussian()), rescor = TRUE)"),
       family = NA, data = d)
})
S("cens() x cosy()", function() {
  d <- eq_dat(seed = 154)
  d$c <- ifelse(d$y > 1.5, "right", "none")
  d$y <- pmin(d$y, 1.5)
  list(bf = "bf(y | cens(c) ~ x + cosy(gr = g))", data = d)
})
S("cumulative x bar_crossing", function() list(
  bf = "bf(y ~ x * (1 | g))", family = "cumulative()",
  data = eq_dat("ordinal", seed = 155)))
S("dirichlet x ar()", function() {
  d <- eq_base(seed = 156)
  n <- nrow(d)
  a <- cbind(2, 2 * exp(0.3 * d$x), 2 * exp(-0.2 + 0.3 * d$z))
  y <- t(apply(a, 1, function(r) {
    v <- rgamma(3, r)
    v / sum(v)
  }))
  colnames(y) <- c("a", "b", "c")
  d$y <- y
  list(bf = "bf(y ~ x + ar(time = t, gr = g))", family = "dirichlet()",
       data = d)
})
S("dirichlet_multinomial x trunc()", function() {
  d <- eq_base(seed = 157)
  n <- nrow(d)
  d$ntr <- 10L
  d$y <- t(vapply(seq_len(n), function(i) {
    stats::rmultinom(1, 10, c(0.3, 0.3, 0.4))[, 1]
  }, numeric(3)))
  colnames(d$y) <- c("a", "b", "c")
  list(bf = "bf(y | trials(ntr) + trunc(ub = 10) ~ x)",
       family = "dirichlet_multinomial()", data = d)
})
S("fcor() x nl", function() {
  d <- eq_dat(seed = 158)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  M <- diag(nrow(d))
  # brms's spelling: in the nl body brms reads fcor(M) as R code and
  # fails on M (review 2026-10-28, finding 1); autocor = is where brms
  # takes a residual correlation term of a nonlinear model
  list(bf = paste0("bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, ",
                   "autocor = ~ fcor(M))"),
       data = d, data2 = list(M = M))
})
S("fcor() x rescor", function() {
  d <- eq_mv(159)
  M <- diag(nrow(d))
  list(bf = "mvbf(bf(y1 ~ x + fcor(M)), bf(y2 ~ x), rescor = TRUE)",
       data = d, data2 = list(M = M))
})
S("gaussian x bar_crossing", function() list(
  bf = "bf(y ~ x * (1 | g))", data = eq_dat(seed = 160)))
S("gaussian x cs_pred()", function() list(
  bf = "bf(y ~ x + cs(z))", data = eq_dat(seed = 161)))
S("gaussian x disc", function() list(
  bf = "bf(y ~ x, disc ~ z)", data = eq_dat(seed = 162)))
S("hurdle_gamma x bar_crossing", function() list(
  bf = "bf(y ~ x * (1 | g))", family = "hurdle_gamma()",
  data = eq_dat("hurdle_gamma", seed = 163)))
S("hurdle_lognormal x fcor()", function() {
  d <- eq_dat("hurdle_lognormal", seed = 164)
  list(bf = "bf(y ~ x + fcor(M))", family = "hurdle_lognormal()", data = d,
       data2 = list(M = diag(nrow(d))))
})
S("hurdle_poisson x bar_crossing", function() list(
  bf = "bf(y ~ x * (1 | g))", family = "hurdle_poisson()",
  data = eq_dat("hurdle_poisson", seed = 165)))
S("nbinom1 x nl", function() {
  d <- eq_dat(seed = 166)
  d$y <- rnbinom(nrow(d), mu = 2 * exp(0.3 * d$x), size = 2)
  list(bf = "bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)",
       family = "nbinom1()", data = d, nobrms = "brms has no nbinom1 family")
})
S("poisson x arma()", function() list(
  # cov = TRUE: brms refuses an MA part under cov = FALSE for this family
  # and builds the latent-residual model under cov = TRUE (review
  # 2026-10-28, finding 1)
  bf = "bf(y ~ x + arma(time = t, gr = g, cov = TRUE))", family = "poisson()",
  data = eq_dat("poisson", seed = 167)))
S("poisson x cs_pred()", function() list(
  bf = "bf(y ~ x + cs(z))", family = "poisson()",
  data = eq_dat("poisson", seed = 168)))
S("von_mises x mi()", function() {
  d <- eq_dat("von_mises", seed = 169)
  d$y[sample(nrow(d), 15)] <- NA
  list(bf = "bf(y | mi() ~ x)", family = "von_mises()", data = d)
})
S("weibull x cs_pred()", function() list(
  bf = "bf(y ~ x + cs(z))", family = "weibull()",
  data = eq_dat("weibull", seed = 170)))
S("zero_inflated_beta x nl", function() {
  d <- eq_dat("zero_inflated_beta", seed = 171)
  list(bf = "bf(y ~ 1 / (1 + exp(-(b1 + b2 * x))), b1 + b2 ~ 1, nl = TRUE)",
       family = "zero_inflated_beta()", data = d)
})
S("zero_inflated_poisson x equidistant", function() list(
  bf = "bf(y ~ x)", family = "zero_inflated_poisson(threshold = 'equidistant')",
  data = eq_dat("zero_inflated_poisson", seed = 172)))
