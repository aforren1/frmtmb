# Log-density identity against brms.
#
# The claim is that frmtmb's objective is the same function of the
# parameters as the Stan program brms generates with flat priors. It is
# checked at a point, never through the two optimizers and never against
# posterior summaries: the estimand differs and any tolerance loose
# enough to pass would also hide a real divergence.
#
# Three checks, all in brms_lp_check() in helper-brms.R:
#
#   A  log_prob(sf, upars, adjust_transform = FALSE) at frmtmb's
#      estimates equals logLik(fit) plus a known constant, to
#      1e-6 * max(1, |logLik|). The only admitted constants are the flat
#      Dirichlet brms keeps on a mo() simplex and on mixture weights.
#      Any other nonzero constant is a finding, not a tolerance.
#   B  grad_log_prob at that same point has max absolute entry below
#      1e-3. This is what catches a wrong parameter map: a mistranslated
#      parameter lands off brms's optimum even when the value happens to
#      look plausible.
#   C  for random-effect models brms has no marginal likelihood, so A
#      and B run on the joint density at frmtmb's outer estimates plus
#      its conditional modes, against minus RTMB's inner objective plus
#      the map's log-Jacobian, and B is asserted on the z block only.
#
# There is no known-divergence list for numeric mismatches. The only
# admissible exemption is a design choice frmtmb states in its own
# source, recorded here with its reason. That list has one entry, row
# 3's `mo(inc) * z`, and it is recorded by asserting the structural
# difference rather than by skipping the row.
#
# Stan compiles here. The whole file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
# Compiled programs are cached under FRMTMB_STAN_CACHE, keyed by the
# hash of the Stan code and the rstan version, so only the first run of
# a shape pays for it. See dev/brms-likelihood-tests.md.

# ---------------------------------------------------------------------
# The translator, tested on its own before any density is compared.
# ---------------------------------------------------------------------

test_that("the translator round-trips through Stan's constraints", {
  skip_unless_brms_fit()

  set.seed(11)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x))
  bform <- brms::bf(y ~ x + z, sigma ~ x)
  fit <- frm(bf(y ~ x + z, sigma ~ x) + gaussian(), data = dd)

  prior <- brms_flat_prior(bform, data = dd, family = gaussian())
  code <- brms::make_stancode(bform, data = dd, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = dd, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)

  # every parameter brms declares got a rule, and no others
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  # brms centers X inside the Stan program, so its Intercept is not the
  # intercept frmtmb reports; the generated quantity is
  #   b_Intercept = Intercept - dot_product(means_X, b)
  fe <- fixef(fit)
  expect_lt(abs(pars[["Intercept"]] -
                  sum(colMeans(sdat$X)[-1] * pars[["b"]]) -
                  fe$mu[["(Intercept)"]]), 1e-10)
  expect_lt(abs(pars[["Intercept_sigma"]] -
                  sum(colMeans(sdat$X_sigma)[-1] * pars[["b_sigma"]]) -
                  fe$sigma[["(Intercept)"]]), 1e-10)
})

test_that("the simplex and group-level rules round-trip", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
  bform <- brms::bf(y ~ mo(inc) + z)
  fit <- frm(bf(y ~ mo(inc) + z) + gaussian(), data = dm)
  prior <- brms_flat_prior(bform, data = dm, family = gaussian())
  code <- brms::make_stancode(bform, data = dm, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = dm, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_par_roundtrip(sf, pars)
  # the simplex is frmtmb's zeta through a softmax with a fixed zero
  expect_lt(abs(sum(pars[["simo_1"]]) - 1), 1e-12)
  expect_length(pars[["simo_1"]], sdat$Jmo[[1]])
  # sigma has no linear predictor, so brms declares it on the natural
  # scale while frmtmb estimates it through the log link
  expect_lt(abs(log(pars[["sigma"]]) - fixef(fit)$sigma[["(Intercept)"]]),
            1e-12)

  data(sleepstudy, package = "lme4")
  bform2 <- brms::bf(Reaction ~ Days + (Days | Subject), sigma ~ Days)
  fit2 <- frm(bf(Reaction ~ Days + (Days | Subject), sigma ~ Days) +
                gaussian(), data = sleepstudy)
  prior2 <- brms_flat_prior(bform2, data = sleepstudy,
                            family = gaussian())
  code2 <- brms::make_stancode(bform2, data = sleepstudy,
                               family = gaussian(), prior = prior2)
  sdat2 <- brms_standata(bform2, data = sleepstudy, family = gaussian(),
                         prior = prior2)
  sf2 <- suppressMessages(rstan::sampling(brms_stan_model(code2),
                                          data = sdat2, chains = 0))
  rtab2 <- brms_ranef_table(bform2, sleepstudy, gaussian(), prior2)
  pars2 <- stan_pars_from_fit(fit2, sdat2, code2, rtab2)
  expect_par_roundtrip(sf2, pars2)

  # sd and L reproduce frmtmb's covariance, and z inverts brms's
  # r = (diag(sd) L z)^T against ranef() rather than the raw b vector,
  # whose order follows frmtmb's Zt and is not brms's
  sigma_hat <- diag(pars2[["sd_1"]]) %*% tcrossprod(pars2[["L_1"]]) %*%
    diag(pars2[["sd_1"]])
  expect_lt(max(abs(sigma_hat - unclass(VarCorr(fit2))[[1]])), 1e-8)
  r <- t(diag(pars2[["sd_1"]]) %*% pars2[["L_1"]] %*% pars2[["z_1"]])
  info <- brms_group_info(rtab2, 1)
  expect_lt(max(abs(r - ranef(fit2)[[1]][info$labels, ])), 1e-8)
  # brms's coefficient order inside a group is Intercept then Days,
  # matching standata's Z_1_1 then Z_1_2. get_prior() sorts its rows
  # alphabetically instead, which is why the order is not read there.
  expect_identical(info$coefs, c("Intercept", "Days"))
})

# ---------------------------------------------------------------------
# The model matrix. One call per row.
# ---------------------------------------------------------------------

test_that("row 1: distributional gaussian, y ~ x + z with sigma ~ x", {
  skip_unless_brms_fit()

  set.seed(11)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x))
  fit <- frm(bf(y ~ x + z, sigma ~ x) + gaussian(), data = dd)
  brms_lp_check(brms::bf(y ~ x + z, sigma ~ x), gaussian(), dd, fit)
})

test_that("row 2: monotonic effects, y ~ mo(inc) + z", {
  skip_unless_brms_fit()

  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
  fit <- frm(bf(y ~ mo(inc) + z) + gaussian(), data = dm)
  # the one prior the empty string cannot remove is the flat Dirichlet
  # on the simplex, whose normalizing constant is lgamma(D)
  brms_lp_check(brms::bf(y ~ mo(inc) + z), gaussian(), dm, fit,
                const = lgamma(3))
})

test_that("row 3: monotonic interaction, y ~ mo(inc):z", {
  skip_unless_brms_fit()

  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
  fit <- frm(bf(y ~ mo(inc):z) + gaussian(), data = dm)
  brms_lp_check(brms::bf(y ~ mo(inc):z), gaussian(), dm, fit,
                const = lgamma(3))
})

test_that("row 3: mo(inc) * z is a different model in the two packages", {
  skip_unless_brms()

  # EXEMPTION, and the only one. brms builds one simplex per special
  # term, frmtmb one per mo() VARIABLE: R/frame.R keys `mo_zetas` on
  # deparse1(mexpr) and says "simplexes are shared per mo() variable".
  # So mo(inc) and mo(inc):z share a simplex here and have their own in
  # brms, frmtmb's model has two fewer free parameters, and no parameter
  # map can turn one into the other. Nothing is skipped and no tolerance
  # is widened: the structural difference is asserted, so this fails
  # loudly if either package changes its mind. See
  # dev/brms-likelihood-tests.md for which side is right.
  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)

  sdat <- brms_standata(brms::bf(y ~ mo(inc) * z), data = dm,
                        family = gaussian())
  expect_identical(as.integer(sdat$Imo), 2L)
  expect_identical(as.integer(sdat$Jmo), c(3L, 3L))

  fit <- frm(bf(y ~ mo(inc) * z) + gaussian(), data = dm)
  zetas <- grep("^zeta", names(fit$estimates), value = TRUE)
  expect_identical(zetas, "zeta1")
  # both monotonic coefficients are present; only the shape is shared
  expect_true(all(c("moinc", "moinc:z") %in% names(fixef(fit)$mu)))
})

test_that("row 5: nonlinear, y ~ a * exp(-b * x) with a + b ~ 1", {
  skip_unless_brms_fit()

  set.seed(7)
  n <- 120
  dn <- data.frame(x = runif(n, 0, 3))
  dn$y <- 2.5 * exp(-0.8 * dn$x) + rnorm(n, 0, 0.15)
  fit <- frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) + gaussian(),
             data = dn)
  # nlpar predictors are not centered, so every column of X_a and X_b
  # belongs to b_a and b_b and there is no Intercept_a to translate
  brms_lp_check(brms::bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE),
                gaussian(), dn, fit)
})

test_that("row 12: ordinal families, cumulative sratio cratio acat", {
  skip_unless_brms_fit()

  # Ordinal is the one shape where the centering sign flips. brms puts
  # no intercept column in X (Kc == K), declares Intercept as the
  # ordered threshold vector, and its generated quantity is
  #   b_Intercept = Intercept + dot_product(means_X, b)
  # because the thresholds are compared to mu rather than added to it.
  #
  # It is also the one shape where frmtmb's own storage differs between
  # families: cumulative and sratio hold (tau_1, log increments) while
  # cratio and acat hold the thresholds themselves. brms_ord_thresholds()
  # carries that distinction, and these four rows are what pin it down.
  set.seed(5)
  n <- 300
  do <- data.frame(x = rnorm(n))
  do$y <- ordered(cut(0.9 * do$x + rlogis(n),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))

  brms_lp_check(brms::bf(y ~ x), brms::cumulative(), do,
                frm(bf(y ~ x) + cumulative(), data = do))
  brms_lp_check(brms::bf(y ~ x), brms::sratio(), do,
                frm(bf(y ~ x) + sratio(), data = do))
  brms_lp_check(brms::bf(y ~ x), brms::cratio(), do,
                frm(bf(y ~ x) + cratio(), data = do))
  brms_lp_check(brms::bf(y ~ x), brms::acat(), do,
                frm(bf(y ~ x) + acat(), data = do))

  # cs(): brms declares bcs as matrix[Kcs, nthres], frmtmb keeps one
  # bcs<j> vector of length nthres per category-specific covariate
  do$z <- rnorm(n)
  brms_lp_check(brms::bf(y ~ x + cs(z)), brms::sratio(), do,
                frm(bf(y ~ x + cs(z)) + sratio(), data = do))
})

test_that("row 13: categorical(y ~ x)", {
  skip_unless_brms_fit()

  set.seed(31)
  n <- 300
  dc <- data.frame(x = rnorm(n))
  p <- cbind(1, exp(0.3 + 0.5 * dc$x), exp(-0.2 + 0.9 * dc$x))
  p <- p / rowSums(p)
  dc$y <- factor(apply(p, 1, function(pr) sample(1:3, 1, prob = pr)))
  fit <- frm(bf(y ~ x) + categorical(), data = dc)
  brms_lp_check(brms::bf(y ~ x), brms::categorical(), dc, fit)
})

test_that("row 14: censoring and known standard errors", {
  skip_unless_brms_fit()

  set.seed(11)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x))
  dd$cens <- sample(c(0L, 1L), n, TRUE, prob = c(0.8, 0.2))
  dd$s <- runif(n, 0.5, 1.5)

  fit_c <- frm(bf(y | cens(cens) ~ x) + gaussian(), data = dd)
  brms_lp_check(brms::bf(y | cens(cens) ~ x), gaussian(), dd, fit_c)

  fit_s <- frm(bf(y | se(s) ~ x) + gaussian(), data = dd)
  brms_lp_check(brms::bf(y | se(s) ~ x), gaussian(), dd, fit_s)
})

test_that("row 15: binomial with trials(n)", {
  skip_unless_brms_fit()

  set.seed(9)
  n <- 200
  db <- data.frame(x = rnorm(n), n = sample(3:12, n, TRUE))
  db$y <- rbinom(n, db$n, plogis(0.3 + 0.7 * db$x))
  fit <- frm(bf(y | trials(n) ~ x) + binomial(), data = db)
  brms_lp_check(brms::bf(y | trials(n) ~ x), binomial(), db, fit)
})

test_that("row 16: zero-inflated poisson with zi ~ x", {
  skip_unless_brms_fit()

  set.seed(13)
  n <- 300
  dz <- data.frame(x = rnorm(n))
  dz$y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.3 * dz$x)), 0L,
                 rpois(n, exp(0.6 + 0.4 * dz$x)))
  fit <- frm(bf(y ~ x, zi ~ x) + zero_inflated_poisson(), data = dz)
  brms_lp_check(brms::bf(y ~ x, zi ~ x),
                brms::zero_inflated_poisson(), dz, fit)
})

test_that("row 20: weights(w)", {
  skip_unless_brms_fit()

  # brms multiplies each observation's log density by its weight rather
  # than normalizing the weights, and so does frmtmb: the constant is
  # zero, not log of any scaling factor.
  set.seed(11)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x))
  dd$w <- runif(n, 0.5, 2)
  fit <- frm(bf(y | weights(w) ~ x) + gaussian(), data = dd)
  brms_lp_check(brms::bf(y | weights(w) ~ x), gaussian(), dd, fit)
})

test_that("row 21: one-predictor fits across the family roster", {
  skip_unless_brms_fit()

  # A and B only. These carry no special terms, so what they exercise is
  # the family's own log density and the natural-scale rule for a dpar
  # with no linear predictor: shape for Gamma and for negbinomial.
  set.seed(17)
  n <- 250
  d <- data.frame(x = rnorm(n))

  d$y <- rpois(n, exp(0.5 + 0.6 * d$x))
  brms_lp_check(brms::bf(y ~ x), poisson(), d,
                frm(bf(y ~ x) + poisson(), data = d))

  d$y <- rgamma(n, shape = 3, rate = 3 / exp(0.7 + 0.4 * d$x))
  brms_lp_check(brms::bf(y ~ x), Gamma(link = "log"), d,
                frm(bf(y ~ x) + Gamma(link = "log"), data = d))

  d$y <- rnbinom(n, size = 2, mu = exp(0.5 + 0.5 * d$x))
  brms_lp_check(brms::bf(y ~ x), brms::negbinomial(), d,
                frm(bf(y ~ x) + negbinomial(), data = d))

  d$y <- rbinom(n, 1, plogis(0.2 + 0.8 * d$x))
  brms_lp_check(brms::bf(y ~ x), brms::bernoulli(), d,
                frm(bf(y ~ x) + bernoulli(), data = d))
})

test_that("check C: sleepstudy joint density with (Days | Subject)", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # This shape replaces the posterior-mean comparison retired from
  # test-brms-agreement.R: the same model, checked as an identity of
  # densities instead of an inequality between two estimands.
  fit <- frm(bf(Reaction ~ Days + (Days | Subject), sigma ~ Days) +
               gaussian(), data = sleepstudy)
  brms_lp_check(brms::bf(Reaction ~ Days + (Days | Subject),
                         sigma ~ Days),
                gaussian(), sleepstudy, fit, joint = TRUE)
})

test_that("row 17: mixture(gaussian, gaussian) with theta1 ~ x", {
  skip_unless_brms_fit()

  # brms identifies the mixture by declaring ordered_Intercept, whose
  # entries ARE the component intercepts Intercept_mu1 and
  # Intercept_mu2. With a predictor on theta1 there is no simplex and
  # so no Dirichlet: brms uses softmax over the linear predictors. The
  # constant the plan admits for mixture weights belongs to the
  # simplex spelling, and for two components it is lgamma(2), which is
  # zero anyway.
  set.seed(37)
  n <- 400
  dx <- data.frame(x = rnorm(n))
  k <- rbinom(n, 1, 0.35)
  dx$y <- ifelse(k == 1, rnorm(n, 3, 1), rnorm(n, -1, 1))
  fit <- frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()),
             data = dx)
  brms_lp_check(brms::bf(y ~ 1, theta1 ~ x),
                brms::mixture(gaussian(), gaussian()), dx, fit)
})

test_that("check C: row 7, (1 | q | g) merged across mu and sigma", {
  skip_unless_brms_fit()
  skip_if_not_installed("MASS")

  # The hardest group-level shape. brms merges the two dpars' effects
  # into one block of M_1 = 2 correlated coefficients, and frmtmb names
  # the columns "y.mu:(Intercept)" and "y.sigma:(Intercept)" because the
  # bare coefficient name now appears once per linear predictor.
  set.seed(43)
  ng <- 30
  nper <- 12
  n <- ng * nper
  gg <- factor(rep(seq_len(ng), each = nper))
  u <- MASS::mvrnorm(ng, c(0, 0), matrix(c(0.6, 0.25, 0.25, 0.3), 2))
  dd <- data.frame(x = rnorm(n), g = gg)
  dd$y <- 1 + 0.8 * dd$x + u[gg, 1] +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x + u[gg, 2]))
  fit <- frm(bf(y ~ x + (1 | q | g), sigma ~ x + (1 | q | g)) +
               gaussian(), data = dd)
  brms_lp_check(brms::bf(y ~ x + (1 | q | g), sigma ~ x + (1 | q | g)),
                gaussian(), dd, fit, joint = TRUE)
})

test_that("check C: row 3's (1 | g) variant, monotonic interaction", {
  skip_unless_brms_fit()

  # The plan's C extension of row 3, in the spelling the two packages
  # agree on. The admitted Dirichlet survives into the joint density
  # unchanged, because it sits on the simplex and not on the effects.
  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300),
                   g = factor(rep(1:20, 15)))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
  fit <- frm(bf(y ~ mo(inc):z + (1 | g)) + gaussian(), data = dm)
  brms_lp_check(brms::bf(y ~ mo(inc):z + (1 | g)), gaussian(), dm, fit,
                joint = TRUE, const = lgamma(3))
})

test_that("check C: row 16's (1 | g) variant, zero-inflated poisson", {
  skip_unless_brms_fit()

  # The plan's C extension of row 16. A single-coefficient group has no
  # correlation matrix, so brms declares sd_1 and z_1 and no L_1, and
  # the log-Jacobian is n_levels * sum(log(sd)) alone.
  set.seed(13)
  n <- 300
  dz <- data.frame(x = rnorm(n), g = factor(rep(1:20, 15)))
  dz$y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.3 * dz$x)), 0L,
                 rpois(n, exp(0.6 + 0.4 * dz$x)))
  fit <- frm(bf(y ~ x + (1 | g), zi ~ x) + zero_inflated_poisson(),
             data = dz)
  brms_lp_check(brms::bf(y ~ x + (1 | g), zi ~ x),
                brms::zero_inflated_poisson(), dz, fit, joint = TRUE)
})
