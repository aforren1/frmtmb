# hurdle_negbinomial() and zero_one_inflated_beta(), brms's two families
# that compose from families frmtmb already had.
#
# Every reference density here is written from stats::dnbinom(),
# stats::pnbinom() and stats::dbeta() after brms 2.23.0's Stan functions
# (chunks/fun_hurdle_negbinomial.stan, fun_zero_one_inflated_beta.stan),
# not from the package code. Tolerances are ratios: to the double
# epsilon for a density, to the measured log-likelihood for a fit, and to
# the measured standard error for a coefficient. dev/fams-validate.R has
# the same comparisons with every number printed.

# P(0) = hu; y > 0: (1 - hu) NB(y) / P(Y > 0), the tail from pnbinom()
ref_hnb <- function(y, mu, shape, hu) {
  ifelse(y == 0, log(hu),
         log1p(-hu) + stats::dnbinom(y, size = shape, mu = mu, log = TRUE) -
           stats::pnbinom(0, size = shape, mu = mu, lower.tail = FALSE,
                          log.p = TRUE))
}

ref_zoib <- function(y, mu, phi, zoi, coi) {
  ifelse(y == 0, log(zoi) + log1p(-coi),
         ifelse(y == 1, log(zoi) + log(coi),
                log1p(-zoi) + stats::dbeta(y, mu * phi, (1 - mu) * phi,
                                           log = TRUE)))
}

# relative difference, floored at one so a density near zero is judged on
# its absolute error
rel_diff <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))
ULPS <- 64 * .Machine$double.eps

# The standard error of every fixed coefficient, in the order
# unlist(fixef_by_dpar()) gives them, each on its link scale, read off
# vcov(full = TRUE) by confint()'s internal names.
coef_se <- function(fit) {
  fe <- fixef_by_dpar(fit)
  nm <- unlist(lapply(names(fe), function(dp) {
    if (dp == "mu") names(fe[[dp]]) else paste0(dp, "_", names(fe[[dp]]))
  }))
  sqrt(diag(vcov(fit, full = TRUE)))[nm]
}

sim_hnb <- function(seed, n = 500, ngrp = 0, shape = 1.3) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  g <- factor(if (ngrp > 0) rep(seq_len(ngrp), length.out = n) else 1)
  u <- if (ngrp > 0) rnorm(ngrp, 0, 0.5)[g] else 0
  mu <- exp(0.6 + 0.4 * x + u)
  hu <- plogis(-0.4 + 0.7 * z)
  p0 <- stats::dnbinom(0, size = shape, mu = mu)
  yp <- pmax(stats::qnbinom(p0 + runif(n) * (1 - p0), size = shape,
                            mu = mu), 1)
  data.frame(y = ifelse(runif(n) < hu, 0, yp), x = x, z = z, g = g)
}

sim_zoib <- function(seed, n = 600) {
  set.seed(seed)
  x <- rnorm(n)
  mu <- plogis(-0.2 + 0.5 * x)
  zoi <- plogis(-1 + 0.6 * x)
  coi <- plogis(0.3 - 0.8 * x)
  y <- ifelse(runif(n) < zoi, as.numeric(runif(n) < coi),
              stats::rbeta(n, mu * 6, (1 - mu) * 6))
  data.frame(y = y, x = x)
}

# ------------------------------------------------------------ constructors

test_that("the constructors carry brms's names, dpars and default links", {
  f <- hurdle_negbinomial()
  expect_identical(f$family, "hurdle_negbinomial")
  expect_identical(f$dpars, c("mu", "shape", "hu"))
  expect_identical(c(f$link, f$link_shape, f$link_hu),
                   c("log", "log", "logit"))
  expect_identical(hurdle_negbinomial(log)$link, "log")
  g <- zero_one_inflated_beta()
  expect_identical(g$family, "zero_one_inflated_beta")
  expect_identical(g$dpars, c("mu", "phi", "zoi", "coi"))
  expect_identical(c(g$link, g$link_phi, g$link_zoi, g$link_coi),
                   c("logit", "log", "logit", "logit"))
  # brms's family objects and names resolve to the same families
  expect_identical(brmsfamily("hu_negbinomial")$family, "hurdle_negbinomial")
  # the shape of a brms family object, without needing brms
  bf_obj <- structure(list(family = "zero_one_inflated_beta",
                           link = "logit", link_coi = "identity"),
                      class = c("brmsfamily", "family"))
  f2 <- frmtmb:::as_frmtmb_family(bf_obj)
  expect_identical(c(f2$family, f2$link_coi),
                   c("zero_one_inflated_beta", "identity"))
})

test_that("links and responses brms refuses are refused by name", {
  expect_error(hurdle_negbinomial("inverse"),
               paste("'inverse' is not a supported link for family",
                     "'hurdle_negbinomial'"))
  expect_error(zero_one_inflated_beta("sqrt"),
               "not a supported link for family 'zero_one_inflated_beta'")
  expect_error(zero_one_inflated_beta(link_coi = "log"),
               "not a supported link for parameter 'coi'")
  expect_error(hurdle_negbinomial(link_hu = "log"),
               "not a supported link for parameter 'hu'")
  d <- data.frame(y = c(0, 1, 2.5, 3), x = 1:4)
  expect_error(frm(bf(y ~ x) + hurdle_negbinomial(), data = d),
               "hurdle_negbinomial: response must be non-negative integers")
  d <- data.frame(y = c(0, 0.3, 1, 1.2), x = 1:4)
  expect_error(frm(bf(y ~ x) + zero_one_inflated_beta(), data = d),
               "zero_one_inflated_beta: response must be in [0, 1]",
               fixed = TRUE)
  # brms bars this one as a mixture component and not the other
  expect_error(mixture(Beta, zero_one_inflated_beta),
               "not allowed in mixture models: zero_one_inflated_beta",
               fixed = TRUE)
  expect_no_error(mixture(poisson, hurdle_negbinomial))
})

test_that("cens() and trunc() are refused: neither family has a CDF", {
  d <- sim_hnb(5, n = 60)
  expect_error(frm(bf(y | trunc(lb = 1) ~ x) + hurdle_negbinomial(),
                   data = d), "need a family with a CDF")
  expect_error(frm(bf(y | cens(x > 1) ~ x) + hurdle_negbinomial(),
                   data = d), "need a family with a CDF")
})

# ----------------------------------------------------------------- density

test_that("hurdle_negbinomial's density is brms's, on and off the tape", {
  fam <- hurdle_negbinomial()
  y <- c(0, 1, 2, 5, 40)
  # eta = -25 is where 1 - P(0) formed by subtraction keeps five digits
  for (eta in c(-25, -8, 0, 4)) {
    base <- list(mu = exp(eta), shape = 0.7, hu = 0.3)
    tape <- c(base, list(.eta_mu = eta, .eta_shape = log(0.7),
                         .eta_hu = stats::qlogis(0.3)))
    ref <- ref_hnb(y, exp(eta), 0.7, 0.3)
    for (dp in list(base, tape)) {
      got <- as.numeric(fam$lpdf(y, lapply(dp, rep, length.out = 5L),
                                 list()))
      expect_lt(rel_diff(got, ref), ULPS)
    }
  }
  # a mean link other than log goes through log(mu) itself
  fs <- hurdle_negbinomial(link = "sqrt")
  got <- as.numeric(fs$lpdf(y, list(mu = rep(2.5, 5), shape = rep(1.5, 5),
                                    hu = rep(0.2, 5)), list()))
  expect_lt(rel_diff(got, ref_hnb(y, 2.5, 1.5, 0.2)), ULPS)
})

test_that("zero_one_inflated_beta's density is brms's, atoms and interior", {
  fam <- zero_one_inflated_beta()
  y <- c(0, 1, 0.2, 0.5, 0.97)
  dp <- list(mu = 0.35, phi = 4, zoi = 0.25, coi = 0.6)
  got <- as.numeric(fam$lpdf(y, lapply(dp, rep, length.out = 5L), list()))
  ref <- ref_zoib(y, 0.35, 4, 0.25, 0.6)
  expect_lt(rel_diff(got, ref), ULPS)
  tape <- c(dp, list(.eta_mu = stats::qlogis(0.35), .eta_phi = log(4),
                     .eta_zoi = stats::qlogis(0.25),
                     .eta_coi = stats::qlogis(0.6)))
  got <- as.numeric(fam$lpdf(y, lapply(tape, rep, length.out = 5L),
                             list()))
  expect_lt(rel_diff(got, ref), ULPS)
})

# ------------------------------------------------------------ fit agreement

test_that("hurdle_negbinomial matches glmmTMB's truncated_nbinom2 hurdle", {
  skip_if_not_installed("glmmTMB")
  d <- sim_hnb(101)
  fit <- frm(bf(y ~ x, hu ~ z) + hurdle_negbinomial(), data = d)
  # the objective at the optimum against the hand-written likelihood
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  ll <- as.numeric(stats::logLik(fit))
  expect_lt(abs(ll - sum(ref_hnb(d$y, dp$mu, dp$shape, dp$hu))),
            ULPS * abs(ll))
  ref <- suppressWarnings(glmmTMB::glmmTMB(
    y ~ x, ziformula = ~z, family = glmmTMB::truncated_nbinom2, data = d))
  # optimizer precision: both sit at the same maximum, so the gap in
  # logLik is second order in the gap in the estimates
  expect_lt(abs(ll - as.numeric(stats::logLik(ref))), 1e-8 * abs(ll))
  fe <- fixef_by_dpar(fit)
  se <- coef_se(fit)
  expect_lt(max(abs(unlist(fe) - c(glmmTMB::fixef(ref)$cond,
                                   log(stats::sigma(ref)),
                                   glmmTMB::fixef(ref)$zi)) / se), 1e-3)
})

test_that("hurdle_negbinomial with a random effect matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  d <- sim_hnb(202, n = 800, ngrp = 25)
  fit <- frm(bf(y ~ x + (1 | g), hu ~ z) + hurdle_negbinomial(), data = d)
  ref <- suppressWarnings(glmmTMB::glmmTMB(
    y ~ x + (1 | g), ziformula = ~z, family = glmmTMB::truncated_nbinom2,
    data = d))
  ll <- as.numeric(stats::logLik(fit))
  expect_lt(abs(ll - as.numeric(stats::logLik(ref))), 1e-8 * abs(ll))
  fe <- fixef_by_dpar(fit)
  expect_lt(max(abs(fe$mu - glmmTMB::fixef(ref)$cond) /
                  coef_se(fit)[1:2]), 1e-3)
})

test_that("zero_one_inflated_beta factorizes into its three parts", {
  skip_if_not_installed("glmmTMB")
  # IDENTITY: with separate predictors the log-likelihood is a bernoulli
  # on 1{y in {0, 1}}, a bernoulli on y over the boundary rows and a beta
  # over the interior rows, sharing no parameter. The joint ML fit is the
  # three separate ML fits and its logLik is their sum; what is checked
  # is that the family implements that algebra, at optimizer precision.
  d <- sim_zoib(303)
  fit <- frm(bf(y ~ x, zoi ~ x, coi ~ x) + zero_one_inflated_beta(),
             data = d)
  d$edge <- as.numeric(d$y == 0 | d$y == 1)
  r_zoi <- stats::glm(edge ~ x, family = stats::binomial, data = d)
  r_coi <- stats::glm(y ~ x, family = stats::binomial,
                      data = d[d$edge == 1, ])
  r_beta <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::beta_family(),
                             data = d[d$edge == 0, ])
  ll <- as.numeric(stats::logLik(fit))
  ll3 <- as.numeric(stats::logLik(r_zoi)) + as.numeric(stats::logLik(r_coi)) +
    as.numeric(stats::logLik(r_beta))
  expect_lt(abs(ll - ll3), 1e-8 * abs(ll))
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  expect_lt(abs(ll - sum(ref_zoib(d$y, dp$mu, dp$phi, dp$zoi, dp$coi))),
            ULPS * abs(ll))
  fe <- fixef_by_dpar(fit)
  ref3 <- c(glmmTMB::fixef(r_beta)$cond, log(stats::sigma(r_beta)),
            stats::coef(r_zoi), stats::coef(r_coi))
  expect_identical(names(fe), c("mu", "phi", "zoi", "coi"))
  expect_lt(max(abs(unlist(fe) - ref3) / coef_se(fit)), 1e-3)
})

# ---------------------------------------------------------------- post-fit

test_that("fitted() is brms's posterior_epred, and dpar = reaches each gate", {
  d <- sim_hnb(7, n = 200)
  fit <- frm(bf(y ~ x, hu ~ z) + hurdle_negbinomial(), data = d)
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  # brms:::posterior_epred_hurdle_negbinomial(), written out
  m <- with(dp, mu / (1 - (shape / (mu + shape))^shape) * (1 - hu))
  expect_lt(rel_diff(fitted(fit)[, "Estimate"], m), ULPS)
  expect_lt(rel_diff(fitted(fit, dpar = "hu")[, "Estimate"], dp$hu), ULPS)
  expect_lt(rel_diff(as.numeric(frm_linpred(fit, type = "zprob")), dp$hu),
            ULPS)

  dz <- sim_zoib(8, n = 300)
  fz <- frm(bf(y ~ x, zoi ~ x) + zero_one_inflated_beta(), data = dz)
  dp <- frmtmb:::eval_dpars(fz)[[1]]
  m <- with(dp, zoi * coi + mu * (1 - zoi))
  expect_lt(rel_diff(fitted(fz)[, "Estimate"], m), ULPS)
  expect_lt(rel_diff(fitted(fz, dpar = "zoi")[, "Estimate"], dp$zoi), ULPS)
  expect_lt(rel_diff(fitted(fz, dpar = "coi")[, "Estimate"],
                     rep(dp$coi, length.out = nrow(dz))), ULPS)
  # two gates, so glmmTMB's one-gate vocabulary is refused by name
  expect_error(frm_linpred(fz, type = "zprob"),
               "has two. Use dpar = \"zoi\"", fixed = TRUE)
})

test_that("the variance functions are the densities' own second moments", {
  # summed and integrated from the reference densities, so a slip in
  # either closed form shows here, and pearson residuals divide by it
  fam <- hurdle_negbinomial()
  dp <- list(mu = 2.2, shape = 0.9, hu = 0.35)
  y <- 0:3000
  p <- exp(ref_hnb(y, 2.2, 0.9, 0.35))
  m1 <- sum(y * p)
  expect_lt(rel_diff(fam$post$mean_fn(dp, list()), m1), 1e-10)
  expect_lt(rel_diff(fam$post$var_fn(dp, list()), sum(y^2 * p) - m1^2),
            1e-10)
  fz <- zero_one_inflated_beta()
  dp <- list(mu = 0.4, phi = 5, zoi = 0.2, coi = 0.7)
  f <- function(k) {
    stats::integrate(function(y) y^k * exp(ref_zoib(y, 0.4, 5, 0.2, 0.7)),
                     0, 1, rel.tol = 1e-12)$value + 0.2 * 0.7
  }
  m1 <- f(1)
  expect_lt(rel_diff(fz$post$mean_fn(dp, list()), m1), 1e-10)
  expect_lt(rel_diff(fz$post$var_fn(dp, list()), f(2) - m1^2), 1e-10)
})

test_that("simulate(), residuals() and conditional_effects() run", {
  d <- sim_hnb(9, n = 200)
  fit <- frm(bf(y ~ x) + hurdle_negbinomial(), data = d)
  s <- simulate(fit, nsim = 2, seed = 1)
  expect_identical(dim(s), c(200L, 2L))
  expect_true(all(s$sim_1 >= 0 & s$sim_1 == round(s$sim_1)))
  expect_identical(nrow(residuals(fit, type = "pearson")), 200L)
  ce <- conditional_effects(fit)
  expect_s3_class(ce, "frmtmb_conditional_effects")

  dz <- sim_zoib(10, n = 200)
  fz <- frm(bf(y ~ x) + zero_one_inflated_beta(), data = dz)
  s <- simulate(fz, nsim = 1, seed = 2)
  expect_true(all(s$sim_1 >= 0 & s$sim_1 <= 1))
  expect_true(any(s$sim_1 == 0) && any(s$sim_1 == 1))
  expect_identical(nrow(residuals(fz, type = "pearson")), 200L)
  ce <- conditional_effects(fz, dpar = "zoi")
  expect_true(all(ce[[1]]$estimate__ > 0 & ce[[1]]$estimate__ < 1))
})

test_that("each new dpar is a prior class, as in brms", {
  d <- sim_zoib(11, n = 150)
  tab <- default_prior(bf(y ~ x, zoi ~ x) + zero_one_inflated_beta(),
                       data = d)
  expect_true(all(c("phi", "coi") %in% tab$class))
  expect_true("zoi" %in% tab$dpar)
  fit <- frm(bf(y ~ x) + zero_one_inflated_beta(), data = d,
             prior = set_prior("beta(2, 2)", class = "coi"))
  expect_identical(prior_summary(fit)$class, "coi")
  dh <- sim_hnb(12, n = 150)
  tab <- default_prior(bf(y ~ x) + hurdle_negbinomial(), data = dh)
  expect_true(all(c("shape", "hu") %in% tab$class))
})

test_that("the compatibility registry lists both families", {
  ft <- frm_compat_features()
  expect_true(all(c("hurdle_negbinomial", "zero_one_inflated_beta") %in%
                    ft$key[ft$kind == "family"]))
  for (f in c("hurdle_negbinomial", "zero_one_inflated_beta")) {
    expect_identical(frm_compat(f, "simulate")$status, "works", info = f)
  }
  expect_identical(frm_compat("hurdle_negbinomial", "cens()")$status,
                   "refused")
})

test_that("the registry says simulate() works where a simulator exists", {
  # PIN: through 0.63.0 the registry's no_simulator group listed three
  # families whose simulators work, so frm_compat() called simulate()
  # refused for them. Seen failing on the base build (dev/fams-findings.md)
  for (f in c("hurdle_poisson", "compois", "tweedie")) {
    expect_identical(frm_compat(f, "simulate")$status, "works", info = f)
    fam <- frmtmb:::family_registry[[f]]()
    expect_true(frmtmb:::sim_can(fam), info = f)
  }
  expect_identical(frm_compat("cox", "simulate")$status, "refused")
})

test_that("one-step-ahead residuals are refused by name on a point mass", {
  # PIN: through 0.63.0 every zero-inflated and hurdle family failed here
  # with base R's "comparison (==) is possible only for atomic and list
  # types". Seen failing on the base build for hurdle_poisson.
  d <- sim_hnb(13, n = 80)
  for (fam in list(hurdle_poisson(), hurdle_negbinomial())) {
    fit <- frm(bf(y ~ x), family = fam, data = d)
    expect_error(residuals(fit, type = "osa"),
                 paste0("not available for family '", fam$family, "'"),
                 fixed = TRUE)
  }
  fz <- frm(bf(y ~ x) + zero_one_inflated_beta(), data = sim_zoib(14, 80))
  expect_error(residuals(fz, type = "osa"), "point mass", fixed = TRUE)
  expect_identical(frm_compat("hurdle_negbinomial", "residuals_osa")$status,
                   "refused")
})

test_that("a response that cannot identify coi warns, and only then", {
  set.seed(15)
  n <- 150
  x <- rnorm(n)
  yb <- stats::rbeta(n, 2, 3)
  edge <- runif(n) < 0.2
  cases <- list(
    "0s but no 1" = ifelse(edge, 0, yb),
    "1s but no 0" = ifelse(edge, 1, yb),
    "no exact 0 or 1" = yb)
  for (k in names(cases)) {
    d <- data.frame(x = x, y = cases[[k]])
    expect_warning(frm(bf(y ~ x) + zero_one_inflated_beta(), data = d),
                   k, fixed = TRUE)
    # a coi held at a constant has nothing to identify
    expect_no_warning(frm(bf(y ~ x, coi = 0.5) + zero_one_inflated_beta(),
                          data = d))
  }
  d <- data.frame(x = x, y = ifelse(edge, rbinom(n, 1, 0.5), yb))
  expect_true(any(d$y == 0) && any(d$y == 1))
  expect_no_warning(frm(bf(y ~ x) + zero_one_inflated_beta(), data = d))
})
