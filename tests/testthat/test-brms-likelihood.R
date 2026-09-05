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
# source, recorded here with its reason, and every one is recorded by
# asserting the structural difference rather than by skipping the row.
# That list held one entry, row 3's `mo(inc) * z`, until every mo()
# TERM was given its own simplex; that row is now an identity. Four
# entries remain, over three rows: the exact `gp()` nugget (row 10a),
# brms's `ar(cov = FALSE)` likelihood (row 18d), and the esicar and
# bym2 CAR parameterizations (row 19c).
#
# Stan compiles here. The whole file is opt-in, and skip_unless_brms()
# calls skip_on_cran(), so outside R CMD check BOTH are needed:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
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

test_that("row 3: mo(inc) * z is the same model in both packages", {
  skip_unless_brms_fit()

  # This row used to hold the file's only exemption. frmtmb keyed the
  # simplex on the mo() VARIABLE, so mo(inc) and mo(inc):z shared one
  # shape, its model had two fewer free parameters than brms's, and no
  # parameter map could join them; the difference was asserted rather
  # than skipped. Every mo() TERM now has its own simplex, so the row
  # is an identity like every other and the exemption list is empty.
  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)

  sdat <- brms_standata(brms::bf(y ~ mo(inc) * z), data = dm,
                        family = gaussian())
  expect_identical(as.integer(sdat$Imo), 2L)
  expect_identical(as.integer(sdat$Jmo), c(3L, 3L))

  fit <- frm(bf(y ~ mo(inc) * z) + gaussian(), data = dm)
  # zeta<j> is simo_<j>: two simplexes, in brms's special-term order
  expect_identical(grep("^zeta", names(fit$estimates), value = TRUE),
                   c("zeta1", "zeta2"))
  expect_true(all(c("moinc", "moinc:z") %in% names(fixef(fit)$mu)))
  # the admitted constant is the flat Dirichlet on EACH simplex
  brms_lp_check(brms::bf(y ~ mo(inc) * z), gaussian(), dm, fit,
                const = 2 * lgamma(3))
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

# ---------------------------------------------------------------------
# Multivariate responses. brms's suffix and frmtmb's prefix are the same
# name in the opposite order: brms declares sigma_y1 where frmtmb keys
# y1_sigma, and family(fit)$links is empty for a multivariate fit
# because each response carries its own family.
# ---------------------------------------------------------------------

# Two correlated responses on 20 groups, with the second built from the
# first so that the residual correlation is not zero and Lrescor is
# actually identified.
brms_mv_data <- function() {
  set.seed(21)
  n <- 120
  dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 6)))
  u <- matrix(rnorm(40, 0, 0.7), 20, 2)
  dd$y1 <- 1 + 0.8 * dd$x + u[dd$g, 1] + rnorm(n, 0, 0.9)
  dd$y2 <- -0.5 + 0.4 * dd$x + u[dd$g, 2] + rnorm(n, 0, 0.6) +
    0.5 * dd$y1
  dd
}

test_that("row 8: mvbf(y1 ~ x, y2 ~ x) with set_rescor(TRUE)", {
  skip_unless_brms_fit()

  dd <- brms_mv_data()
  bform <- brms::mvbf(brms::bf(y1 ~ x), brms::bf(y2 ~ x)) +
    brms::set_rescor(TRUE)
  fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE), data = dd,
             family = gaussian())

  prior <- brms_flat_prior(bform, data = dd, family = gaussian())
  code <- brms::make_stancode(bform, data = dd, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = dd, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  # the two spellings of one name: brms's b_y1 is frmtmb's y1_mu, and
  # brms's sigma_y1 is frmtmb's y1_sigma on the natural scale
  fe <- fixef(fit)
  expect_lt(abs(pars[["b_y1"]][[1]] - fe$y1_mu[["x"]]), 1e-12)
  expect_lt(abs(log(pars[["sigma_y1"]]) - fe$y1_sigma[["(Intercept)"]]),
            1e-12)
  expect_lt(abs(log(pars[["sigma_y2"]]) - fe$y2_sigma[["(Intercept)"]]),
            1e-12)
  # Lrescor is the Cholesky factor of frmtmb's residual correlation, in
  # brms's response order rather than alphabetically
  expect_identical(brms_resp_order(sdat), c("y1", "y2"))
  expect_lt(max(abs(tcrossprod(pars[["Lrescor"]]) -
                      unname(rescor_matrix(fit)))), 1e-10)
  expect_true(all(pars[["Lrescor"]][upper.tri(pars[["Lrescor"]])] == 0))

  brms_lp_check(bform, gaussian(), dd, fit)
})

test_that("check C: row 9, multivariate with (1 | p | g) in both", {
  skip_unless_brms_fit()

  # The group-level block spans two RESPONSES rather than two dpars, so
  # frmtmb repeats the coefficient name once per response and names the
  # columns "y1.mu:(Intercept)" and "y2.mu:(Intercept)". Row 7's dpar
  # fallback cannot separate those two, which is why brms_frm_coef()
  # takes the response as well.
  dd <- brms_mv_data()
  bform <- brms::mvbf(brms::bf(y1 ~ x + (1 | p | g)),
                      brms::bf(y2 ~ x + (1 | p | g))) +
    brms::set_rescor(TRUE)
  fit <- frm(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
               set_rescor(TRUE), data = dd, family = gaussian())

  prior <- brms_flat_prior(bform, data = dd, family = gaussian())
  code <- brms::make_stancode(bform, data = dd, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = dd, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  rtab <- brms_ranef_table(bform, dd, gaussian(), prior)
  pars <- stan_pars_from_fit(fit, sdat, code, rtab)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  # the merged block is one brms id carrying both responses, and the
  # ranef table is where the response of each row comes from
  info <- brms_group_info(rtab, 1)
  expect_identical(info$coefs, c("Intercept", "Intercept"))
  expect_identical(info$resps, c("y1", "y2"))
  expect_identical(as.integer(sdat$M_1), 2L)
  sigma_hat <- diag(pars[["sd_1"]]) %*% tcrossprod(pars[["L_1"]]) %*%
    diag(pars[["sd_1"]])
  expect_lt(max(abs(sigma_hat - unclass(VarCorr(fit))[[1]])), 1e-8)

  brms_lp_check(bform, gaussian(), dd, fit, joint = TRUE)
})

test_that("check C: row 4, mi() imputation and mi(sdx) measurement error", {
  skip_unless_brms_fit()

  # Two latent blocks in one model: y is missing for 30 rows, which brms
  # declares as Ymi_y, and x is observed with a known measurement SD, so
  # every one of its 200 values is latent and brms declares Yl_x. frmtmb
  # keeps both in one `miss` vector and integrates it out, which makes
  # this a joint-density row: logLik(fit) is already marginal.
  set.seed(71)
  n <- 200
  z <- rnorm(n)
  x <- rnorm(n, 0.5 + 0.8 * z, 0.7)
  y <- rnorm(n, 1 + 0.6 * x + 0.3 * z, 0.9)
  ym <- y
  ym[sort(sample(n, 30))] <- NA
  d4 <- data.frame(y = ym, x = x, z = z, sdx = runif(n, 0.2, 0.5))

  bform <- brms::bf(y | mi() ~ mi(x) + z) + brms::bf(x | mi(sdx) ~ z) +
    brms::set_rescor(FALSE)
  fit <- frm(bf(y | mi() ~ mi(x) + z) + gaussian() +
               bf(x | mi(sdx) ~ z) + gaussian(), data = d4)

  prior <- brms_flat_prior(bform, data = d4, family = gaussian())
  code <- brms::make_stancode(bform, data = d4, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = d4, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  # the imputed block covers exactly the missing rows, the latent block
  # covers every row, and the two together are frmtmb's `miss`
  expect_identical(as.integer(sdat$Jmi_y), which(is.na(d4$y)))
  expect_length(pars[["Ymi_y"]], 30L)
  expect_length(pars[["Yl_x"]], n)
  expect_length(fit$estimates[["miss"]], 30L + n)
  expect_lt(max(abs(c(pars[["Ymi_y"]], pars[["Yl_x"]]) -
                      fit$estimates[["miss"]])), 1e-12)
  # mi(x) is a special term, so its coefficient is bsp_y and not a
  # column of X_y
  expect_lt(abs(pars[["bsp_y"]][[1]] - fixef(fit)$y_mu[["mix"]]), 1e-12)

  brms_lp_check(bform, gaussian(), d4, fit, joint = TRUE)
})

test_that("check C: row 6, nonlinear with a ~ 1 + (1 | g)", {
  skip_unless_brms_fit()

  # Nothing in the map is new here: the nlpar rules come from row 5 and
  # the group-level rules from the C anchors. What this row adds is the
  # two together, where brms names the coefficients b_a and b_b and puts
  # the group-level block under the nlpar rather than under a dpar.
  set.seed(7)
  n <- 180
  d6 <- data.frame(x = runif(n, 0, 3), g = factor(rep(1:20, 9)))
  ug <- rnorm(20, 0, 0.4)
  d6$y <- (2.5 + ug[d6$g]) * exp(-0.8 * d6$x) + rnorm(n, 0, 0.15)

  bform <- brms::bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1,
                    nl = TRUE)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) +
               gaussian(), data = d6)

  prior <- brms_flat_prior(bform, data = d6, family = gaussian())
  code <- brms::make_stancode(bform, data = d6, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = d6, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  rtab <- brms_ranef_table(bform, d6, gaussian(), prior)
  pars <- stan_pars_from_fit(fit, sdat, code, rtab)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  # the block belongs to nlpar a, and frmtmb labels it "a: 1 | g"
  info <- brms_group_info(rtab, 1)
  expect_identical(info$nlpars, "a")
  expect_identical(names(ranef(fit)), "a: 1 | g")
  # nlpar predictors are not centered, so b_a carries the intercept
  # itself rather than a centered one
  expect_lt(abs(pars[["b_a"]][[1]] - fixef(fit)$a[["(Intercept)"]]),
            1e-12)

  brms_lp_check(bform, gaussian(), d6, fit, joint = TRUE)
})

# ---------------------------------------------------------------------
# The wiggly rows. brms standardizes every latent block it declares and
# multiplies it up inside the program, so a smooth and a GP arrive here
# as check C shapes just as a group-level effect does.
# ---------------------------------------------------------------------

brms_smooth_data <- function() {
  set.seed(19)
  n <- 200
  d <- data.frame(x = runif(n, -3, 3), z = runif(n, -2, 2))
  d$y <- sin(d$x) + 0.5 * d$z^2 + 0.4 * d$x * d$z + rnorm(n, 0, 0.4)
  d
}

test_that("check C: row 11, s(x) is a random effect with one basis", {
  skip_unless_brms_fit()

  d11 <- brms_smooth_data()
  bform <- brms::bf(y ~ s(x))
  fit <- frm(bf(y ~ s(x)) + gaussian(), data = d11)

  prior <- brms_flat_prior(bform, data = d11, family = gaussian())
  code <- brms::make_stancode(bform, data = d11, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = d11, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  # The two s() bases are NOT the same columns. brms calls
  # mgcv::smoothCon() with diagonal.penalty = TRUE and frmtmb does not,
  # which test-brms-agreement.R already records as a convention
  # divergence: same span, different rotation. Measured here, that
  # rotation is a PERMUTATION, so it is orthogonal, the i.i.d. prior on
  # the coefficients survives it, and the two models are one model
  # written in two orders. A non-orthogonal map would not be, and
  # brms_basis_map() refuses one.
  expect_identical(as.integer(sdat$nb_1), 1L)
  expect_length(pars[["zs_1_1"]], sdat$knots_1[[1]])
  bk <- brms_smooth_term(fit, "mu", 1)[[1]]
  amat <- brms_basis_map(sdat$Zs_1_1, brms_smooth_z(fit, "mu", bk))
  expect_lt(max(abs(crossprod(amat) - diag(ncol(amat)))), 1e-8)
  expect_lt(max(abs(amat - round(amat))), 1e-10)
  expect_lt(max(abs(rowSums(abs(amat)) - 1)), 1e-10)
  # brms builds s = sds * zs from the permuted coefficients
  expect_lt(max(abs(pars[["sds_1"]][[1]] * pars[["zs_1_1"]] -
                      as.numeric(amat %*% brms_block_b(fit, bk)))), 1e-10)
  expect_lt(abs(pars[["sds_1"]][[1]]^2 -
                  unclass(VarCorr(fit))[["s(x)"]][1, 1]), 1e-10)

  # the unpenalized column is one column on both sides, in the same
  # direction on a different scale, so `bs` is frmtmb's coefficient
  # times that ratio. It carries no prior, so the scale is free.
  expect_identical(as.integer(sdat$Ks), 1L)
  xb <- as.numeric(sdat$Xs)
  xf <- as.matrix(brms_lp_of(fit, "mu")$X)[, "s(x).fx1"]
  expect_lt(max(abs(xf / xb - xf[[1]] / xb[[1]])), 1e-10)
  expect_lt(abs(pars[["bs"]][[1]] - fixef(fit)$mu[["s(x).fx1"]] *
                  (xf[[1]] / xb[[1]])), 1e-10)

  brms_lp_check(bform, gaussian(), d11, fit, joint = TRUE)
})

test_that("check C: row 11, t2(x, z) is three bases under one term", {
  skip_unless_brms_fit()

  # The tensor product is where the per-basis indexing earns its keep:
  # brms declares one sds_1 of length nb_1 = 3 and three zs_1_<j>, while
  # frmtmb makes each basis its own block and gives all three the same
  # term label. The blocks are matched by position within the term.
  d11 <- brms_smooth_data()
  bform <- brms::bf(y ~ t2(x, z))
  fit <- frm(bf(y ~ t2(x, z)) + gaussian(), data = d11)

  prior <- brms_flat_prior(bform, data = d11, family = gaussian())
  code <- brms::make_stancode(bform, data = d11, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = d11, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  expect_identical(as.integer(sdat$nb_1), 3L)
  expect_length(pars[["sds_1"]], 3L)
  bks <- brms_smooth_term(fit, "mu", 1)
  expect_length(bks, 3L)
  expect_identical(vapply(bks, `[[`, 0, "dim"),
                   as.numeric(sdat$knots_1))
  for (j in 1:3) {
    expect_lt(max(abs(pars[["sds_1"]][[j]] * pars[[paste0("zs_1_", j)]] -
                        brms_block_b(fit, bks[[j]]))), 1e-10)
  }

  brms_lp_check(bform, gaussian(), d11, fit, joint = TRUE)
})

brms_gp_data <- function() {
  set.seed(23)
  n <- 80
  d <- data.frame(x = sort(runif(n, 0, 10)))
  d$y <- sin(d$x) + rnorm(n, 0, 0.3)
  d
}

test_that("check C: row 10, gp(x, k = 10) Hilbert-space approximation", {
  skip_unless_brms_fit()

  d10 <- brms_gp_data()
  bform <- brms::bf(y ~ gp(x, k = 10))
  fit <- frm(bf(y ~ gp(x, k = 10)) + gaussian(), data = d10)

  prior <- brms_flat_prior(bform, data = d10, family = gaussian())
  code <- brms::make_stancode(bform, data = d10, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = d10, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)

  # The HSGP basis is fixed data and only the coefficient SDs depend on
  # the parameters, so the map is the diagonal one: zgp is the
  # coefficient vector over the square root of brms's spectral density.
  # Both packages rescale the inputs to unit maximum distance here, so
  # lscale needs no correction, unlike the exact form below.
  bk <- brms_gp_block(fit, "mu", 1)
  expect_identical(bk[["covstruct"]], "hsgp")
  th <- fit$estimates$theta[bk$theta_idx]
  expect_lt(abs(pars[["sdgp_1"]][[1]] - exp(th[[1]])), 1e-12)
  expect_lt(abs(pars[["lscale_1"]][[1]] - exp(th[[2]])), 1e-12)
  spd <- exp(2 * th[[1]]) * sqrt(2 * pi) * exp(th[[2]]) *
    exp(-0.5 * exp(2 * th[[2]]) * as.numeric(sdat$slambda_1))
  expect_lt(max(abs(sqrt(spd) * pars[["zgp_1"]] -
                      brms_block_b(fit, bk))), 1e-8)

  brms_lp_check(bform, gaussian(), d10, fit, joint = TRUE)
})

test_that("row 10: the exact gp() nugget is not brms's", {
  skip_unless_brms()
  skip_if_not_installed("mvtnorm")

  # EXEMPTION. Both packages stabilize the exact GP by adding a nugget
  # to the diagonal, and the two nuggets differ by six orders of
  # magnitude: brms adds an absolute 1e-12 (its gp_exp_quad() writes
  # "cov[n, n] += 1e-12") while frmtmb's gp_corr() adds 1e-6 to the
  # CORRELATION, which is 1e-6 * sdgp^2 on the covariance. Everything
  # else about the kernel agrees to floating point once the length-scale
  # is put on brms's scale, so the joint densities differ by exactly the
  # nugget, and on a kernel whose eigenvalues fall below both of them
  # that difference is hundreds of nats. The marginal likelihood barely
  # moves, which is why the fits agree and the joint densities do not.
  #
  # Nothing is skipped and no tolerance is widened: the two constants
  # are asserted from their two sources, the kernels are compared, and
  # the resulting gap is measured. See dev/brms-likelihood-tests.md.
  d10 <- brms_gp_data()
  fit <- frm(bf(y ~ gp(x)) + gaussian(), data = d10)
  bform <- brms::bf(y ~ gp(x))
  prior <- brms_flat_prior(bform, data = d10, family = gaussian())
  code <- brms::make_stancode(bform, data = d10, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = d10, family = gaussian(),
                        prior = prior)

  # brms's nugget, read off the program it generates
  expect_match(code, "cov[n, n] += 1e-12", fixed = TRUE)

  bk <- brms_gp_block(fit, "mu", 1)
  expect_identical(bk[["covstruct"]], "gp")
  th <- fit$estimates$theta[bk$theta_idx]
  sdgp <- exp(th[[1]])

  # brms scales its GP coordinates to unit maximum distance and frmtmb
  # keeps the data scale, so frmtmb's length-scale is dmax times brms's.
  # That part IS translatable and the rule does it.
  xg <- matrix(as.numeric(sdat$Xgp_1), ncol = as.integer(sdat$Dgp_1))
  expect_lt(abs(max(outer(xg[, 1], xg[, 1], "-")^2) - 1), 1e-10)
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_lt(abs(pars[["lscale_1"]][[1]] - exp(th[[2]]) / sdat$dmax_1),
            1e-12)

  kb <- sdgp^2 * exp(-outer(xg[, 1], xg[, 1], "-")^2 /
                       (2 * pars[["lscale_1"]][[1]]^2))
  diag(kb) <- diag(kb) + 1e-12
  kf <- brms_block_cov(fit, bk)
  # off the diagonal the two kernels are the same matrix
  expect_lt(max(abs((kb - kf)[upper.tri(kb)])), 1e-12)
  # on it they differ by exactly the two nuggets and nothing else
  expect_lt(max(abs(diag(kf - kb) - (1e-6 * sdgp^2 - 1e-12))), 1e-16)

  # and that is enough to move the joint density of the same field by
  # far more than any tolerance this file uses, because the kernel has
  # many eigenvalues below the larger nugget and none below the smaller
  fv <- brms_block_b(fit, bk)
  gap <- mvtnorm::dmvnorm(fv, sigma = kb, log = TRUE) -
    mvtnorm::dmvnorm(fv, sigma = kf, log = TRUE)
  expect_gt(abs(gap), 1)
  expect_gt(sum(eigen(kf / sdgp^2, only.values = TRUE)$values < 1e-6),
            10)
})

# ---------------------------------------------------------------------
# Residual autocorrelation. brms declares these on their natural scales
# and frmtmb on unconstrained ones, and frmtmb's own autocor_natural()
# already spells the result the way brms names it.
# ---------------------------------------------------------------------

brms_ac_data <- function() {
  set.seed(29)
  ng <- 25
  nt <- 8
  n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                  time = rep(seq_len(nt), ng))
  d$x <- rnorm(n)
  e <- unlist(lapply(seq_len(ng), function(i) {
    as.numeric(arima.sim(list(ar = 0.6), nt, sd = 0.6))
  }))
  d$y <- 1 + 0.7 * d$x + e
  d
}

test_that("row 18: ar(p = 1), cosy and unstr residual correlation", {
  skip_unless_brms_fit()

  d18 <- brms_ac_data()

  # ar(): the cov = TRUE spelling, which is the one both packages
  # implement. brms's cov = FALSE default is a different likelihood and
  # is the exemption asserted in the next test.
  fit_ar <- frm(bf(y ~ x + ar(time, gr = g, p = 1, cov = TRUE)) +
                  gaussian(), data = d18)
  bform <- brms::bf(y ~ x + ar(time, gr = g, p = 1, cov = TRUE))
  prior <- brms_flat_prior(bform, data = d18, family = gaussian())
  code <- brms::make_stancode(bform, data = d18, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = d18, family = gaussian(),
                        prior = prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit_ar, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)
  ac <- brms_autocor_of(fit_ar)
  expect_identical(ac$struct, "ar")
  arv <- autocor_natural(fit_ar$estimates$thetaac[ac$theta_idx],
                         ac)[["ar[1]"]]
  expect_lt(abs(pars[["ar"]][[1]] - arv), 1e-14)

  # brms's sigma here is NOT frmtmb's sigma. Its cholesky_cor_ar1()
  # returns the Cholesky factor of T / (1 - ar^2) rather than of the
  # correlation matrix T, so what multiplies it is the INNOVATION SD
  # while frmtmb's sigma is the marginal residual SD. Same model, two
  # scales; without the correction check A misses by 5.28 nats on this
  # design and the gradient is 58.
  expect_match(code, "return cholesky_decompose(mat ./ (1 - ar^2));",
               fixed = TRUE)
  expect_lt(abs(pars[["sigma"]] -
                  exp(fixef(fit_ar)$sigma[["(Intercept)"]]) *
                    sqrt(1 - arv^2)), 1e-12)
  brms_lp_check(bform, gaussian(), d18, fit_ar)

  # cosy(): brms bounds the compound-symmetry correlation on [0, 1] and
  # frmtmb on (-1/(d - 1), 1), so a negative estimate would have no brms
  # counterpart at all. This design is positively correlated.
  fit_cs <- frm(bf(y ~ x + cosy(time, gr = g)) + gaussian(), data = d18)
  expect_gt(autocor_natural(fit_cs$estimates$thetaac,
                            brms_autocor_of(fit_cs))[["cosy"]], 0)
  brms_lp_check(brms::bf(y ~ x + cosy(time, gr = g)), gaussian(), d18,
                fit_cs)

  # unstr(): brms's Lcortime is a Cholesky correlation factor and
  # frmtmb's parameterization is that factor, row-normalized, so the map
  # is the factor itself rather than a recomputed one.
  fit_un <- frm(bf(y ~ x + unstr(time, gr = g)) + gaussian(), data = d18)
  bform_un <- brms::bf(y ~ x + unstr(time, gr = g))
  prior_un <- brms_flat_prior(bform_un, data = d18, family = gaussian())
  code_un <- brms::make_stancode(bform_un, data = d18,
                                 family = gaussian(), prior = prior_un)
  sdat_un <- brms_standata(bform_un, data = d18, family = gaussian(),
                           prior = prior_un)
  sf_un <- suppressMessages(rstan::sampling(brms_stan_model(code_un),
                                            data = sdat_un, chains = 0))
  pars_un <- stan_pars_from_fit(fit_un, sdat_un, code_un)
  expect_par_roundtrip(sf_un, pars_un)
  expect_identical(dim(pars_un[["Lcortime"]]),
                   c(as.integer(sdat_un$n_unique_t),
                     as.integer(sdat_un$n_unique_t)))
  expect_lt(max(abs(rowSums(pars_un[["Lcortime"]]^2) - 1)), 1e-12)
  brms_lp_check(bform_un, gaussian(), d18, fit_un)
})

test_that("row 18: brms's ar(cov = FALSE) is another likelihood", {
  skip_unless_brms()

  # EXEMPTION. brms's ar() defaults to cov = FALSE, the residual
  # REGRESSION form, which conditions on the first observations of each
  # group instead of giving them their stationary distribution. frmtmb
  # implements only the marginal residual-covariance form and refuses
  # the other one by name rather than fitting something else under it.
  # The two are different likelihoods on the same data, so there is no
  # parameter map between them and the row is run on cov = TRUE above.
  d18 <- brms_ac_data()
  expect_error(frm(bf(y ~ x + ar(time, gr = g, p = 1)) + gaussian(),
                   data = d18), "cov = TRUE")
  # brms's two spellings really are two programs: the default declares
  # no correlation factor at all and drops the first observation of
  # each group from the AR recursion
  p0 <- brms_flat_prior(brms::bf(y ~ x + ar(time, gr = g, p = 1)),
                        data = d18, family = gaussian())
  code0 <- brms::make_stancode(brms::bf(y ~ x + ar(time, gr = g, p = 1)),
                               data = d18, family = gaussian(),
                               prior = p0)
  expect_false(grepl("Lcortime", code0, fixed = TRUE))
  expect_true(grepl("J_lag", code0, fixed = TRUE))
})

# ---------------------------------------------------------------------
# CAR fields. brms writes both of its CAR densities unnormalized in the
# field and frmtmb keeps a proper density, so the two differ by a
# closed form in the data alone, computed here from brms's OWN standata.
# ---------------------------------------------------------------------

brms_car_data <- function() {
  w <- matrix(0, 16, 16)
  g <- expand.grid(r = 1:4, c = 1:4)
  for (i in 1:16) {
    for (j in 1:16) {
      if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) w[i, j] <- 1
    }
  }
  dimnames(w) <- list(paste0("L", 1:16), paste0("L", 1:16))
  set.seed(42)
  kmat <- diag(rowSums(w)) - w + matrix(1 / (1e-3 * 16)^2, 16, 16)
  phi <- 1.2 * drop(crossprod(chol(solve(kmat)), rnorm(16)))
  loc <- factor(rep(rownames(w), each = 6), levels = rownames(w))
  d <- data.frame(loc = loc, x = rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] +
    rnorm(nrow(d), 0, 0.5)
  list(d = d, W = w)
}

test_that("check C: row 19, car(escar) and car(icar)", {
  skip_unless_brms_fit()

  s <- brms_car_data()
  d19 <- s$d
  wmat <- s$W

  # escar, the proper CAR. brms declares the field itself, so the map is
  # the identity and the only difference is the normalizer: brms drops
  # -Nloc/2 log(2 pi) and 0.5 log det D from a density that is otherwise
  # frmtmb's, term for term.
  bform <- brms::bf(y ~ x + car(W, gr = loc, type = "escar"))
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "escar")) + gaussian(),
             data = d19, data2 = list(W = wmat))
  prior <- brms_flat_prior(bform, data = d19, family = gaussian(),
                           data2 = list(W = wmat))
  code <- brms::make_stancode(bform, data = d19, family = gaussian(),
                              prior = prior, data2 = list(W = wmat))
  sdat <- brms_standata(bform, data = d19, family = gaussian(),
                        prior = prior, data2 = list(W = wmat))
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)
  bk <- brms_car_block(fit)
  expect_lt(abs(pars[["sdcar"]] - exp(fit$estimates$theta[[1]])), 1e-14)
  expect_lt(abs(pars[["car"]] - car_rho(fit$estimates$theta[[2]])), 1e-14)
  expect_lt(max(abs(pars[["rcar"]] - brms_block_b(fit, bk))), 1e-14)
  # brms's Nneigh is the degree vector, which is what its dropped
  # normalizer is built from
  expect_identical(as.numeric(sdat$Nneigh), as.numeric(bk$aux_car$deg))
  brms_lp_check(bform, gaussian(), d19, fit, joint = TRUE,
                const = brms_car_const(sdat, "escar"),
                data2 = list(W = wmat))

  # icar, the intrinsic one. brms standardizes the field by sdcar and
  # adds a soft sum-to-zero term whose precision is exactly the rank-one
  # term frmtmb folds into its precision matrix, so the quadratic forms
  # agree and only the normalizer differs.
  bform2 <- brms::bf(y ~ x + car(W, gr = loc, type = "icar"))
  fit2 <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
              data = d19, data2 = list(W = wmat))
  prior2 <- brms_flat_prior(bform2, data = d19, family = gaussian(),
                            data2 = list(W = wmat))
  code2 <- brms::make_stancode(bform2, data = d19, family = gaussian(),
                               prior = prior2, data2 = list(W = wmat))
  sdat2 <- brms_standata(bform2, data = d19, family = gaussian(),
                         prior = prior2, data2 = list(W = wmat))
  sf2 <- suppressMessages(rstan::sampling(brms_stan_model(code2),
                                          data = sdat2, chains = 0))
  pars2 <- stan_pars_from_fit(fit2, sdat2, code2)
  expect_setequal(names(pars2), brms_stan_par_names(code2))
  expect_par_roundtrip(sf2, pars2)
  # brms's soft constraint and frmtmb's rank-one term are the same
  # number: normal(sum(zcar) | 0, 0.001 * Nloc) against 1 / (con_sd n)^2
  expect_match(code2, "normal_lpdf(sum(zcar) | 0, 0.001 * Nloc)",
               fixed = TRUE)
  expect_identical(brms_car_block(fit2)$aux_car$con_sd, 0.001)
  brms_lp_check(bform2, gaussian(), d19, fit2, joint = TRUE,
                const = brms_car_const(sdat2, "icar"),
                data2 = list(W = wmat))
})

test_that("row 19: esicar and bym2 have different latent variables", {
  skip_unless_brms()

  # EXEMPTION, two of them, and neither is arithmetic.
  #
  # esicar: brms imposes the sum-to-zero constraint HARD, declaring
  # Nloc - 1 free values and setting the last to minus their sum, and
  # normalizes by (Nloc - 1) log tau. frmtmb imposes it softly, keeps
  # all Nloc values, and normalizes by Nloc log tau plus the log
  # determinant of the constrained precision. The gap between the two
  # densities therefore moves with sdcar and is not a constant, so no
  # parameter map closes it. frmtmb's esicar is in fact its icar: the
  # two fits agree to the last digit.
  #
  # bym2: brms keeps the spatial and the non-spatial parts as SEPARATE
  # latent vectors, 2 * Nloc of them plus rhocar, and frmtmb integrates
  # the mixture into one dense marginal covariance over Nloc values.
  # The two joint densities are functions of different arguments.
  s <- brms_car_data()
  d19 <- s$d
  wmat <- s$W

  for (ty in c("esicar", "bym2")) {
    bform <- brms::bf(stats::as.formula(
      paste0("y ~ x + car(W, gr = loc, type = \"", ty, "\")")))
    prior <- brms_flat_prior(bform, data = d19, family = gaussian(),
                             data2 = list(W = wmat))
    code <- brms::make_stancode(bform, data = d19, family = gaussian(),
                                prior = prior, data2 = list(W = wmat))
    nms <- brms_stan_par_names(code)
    fit <- frm(bf(stats::as.formula(
      paste0("y ~ x + car(W, gr = loc, type = \"", ty, "\")"))) +
        gaussian(), data = d19, data2 = list(W = wmat))
    bk <- brms_car_block(fit)
    expect_identical(bk$aux_car$type, ty)
    expect_length(bk$b_idx, 16L)
    # the translator refuses by name rather than producing a map that
    # cannot exist
    expect_error(stan_pars_from_fit(fit, brms_standata(
      bform, data = d19, family = gaussian(), prior = prior,
      data2 = list(W = wmat)), code), "different set of latent")
    if (identical(ty, "esicar")) {
      # Nloc - 1 free values against frmtmb's Nloc
      expect_true("zcar" %in% nms)
      expect_true(grepl("rcar[Nloc] = - sum(zcar)", code, fixed = TRUE))
    } else {
      # two latent vectors and a mixing proportion against one vector
      expect_true(all(c("zcar", "nszcar", "rhocar") %in% nms))
    }
  }

  # frmtmb's esicar and its icar are the same model, which is the other
  # half of the first divergence: brms's esicar is the constrained one
  # and its icar is not.
  f_ic <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
              data = d19, data2 = list(W = wmat))
  f_es <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
              data = d19, data2 = list(W = wmat))
  expect_lt(abs(as.numeric(logLik(f_ic)) - as.numeric(logLik(f_es))),
            1e-10)
})
