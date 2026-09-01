# Simulation-workflow ergonomics: natural-scale newparams, prior-draw
# (prior-predictive) simulation, censored simulation, and
# conditional_effects() prediction intervals that respect the addition
# terms.

test_that("natural and internal newparams spell the same simulation", {
  set.seed(19)
  dd <- data.frame(x = rnorm(300), g = factor(rep(1:15, 20)), y = 0)
  form <- bf(y ~ x + (1 | g)) + gaussian()

  s_int <- frm_simulate(form, dd, nsim = 3, seed = 1,
                        newparams = list(beta = c(1, 0.5),
                                         betad = log(0.5),
                                         theta = log(0.8)))
  s_nat <- frm_simulate(form, dd, nsim = 3, seed = 1,
                        newparams = list(Intercept = 1, x = 0.5,
                                         sigma = 0.5,
                                         sd_g__Intercept = 0.8))
  expect_equal(s_nat, s_int)
})

test_that("natural newparams round-trip through an us block", {
  set.seed(3)
  dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)), y = 0)
  form <- bf(y ~ x + (1 + x | g)) + gaussian()
  fr <- frm(form, dd, dry_run = "frame")
  slots <- nat_slots(fr)
  expect_true(all(c("Intercept", "x", "sigma", "sd_g__Intercept",
                    "sd_g__x", "cor_g__Intercept__x") %in% names(slots)))

  np <- list(Intercept = 1, x = 0.5, sigma = 0.6,
             sd_g__Intercept = 0.8, sd_g__x = 0.4,
             cor_g__Intercept__x = -0.3)
  est <- apply_natural(fr$par_template, fr, slots, np)
  V <- covstruct_registry$us$vcov(est$theta, fr$re_blocks[[1L]])
  expect_equal(unname(sqrt(diag(V))), c(0.8, 0.4), tolerance = 1e-10)
  expect_equal(unname(stats::cov2cor(V)[1, 2]), -0.3, tolerance = 1e-10)
  expect_equal(unname(est$betad[1]), log(0.6), tolerance = 1e-12)

  # the two spellings still agree once the theta segment is written out
  s_nat <- frm_simulate(form, dd, newparams = np, nsim = 2, seed = 5)
  s_int <- frm_simulate(form, dd, nsim = 2, seed = 5,
                        newparams = list(beta = est$beta,
                                         betad = est$betad,
                                         theta = est$theta))
  expect_equal(s_nat, s_int)
})

test_that("a refit recovers natural-scale simulation parameters", {
  skip_on_cran()
  set.seed(21)
  dd <- data.frame(x = rnorm(600), g = factor(rep(1:40, 15)), y = 0)
  form <- bf(y ~ x + (1 | g)) + gaussian()
  s <- frm_simulate(form, dd, nsim = 1, seed = 4,
                    newparams = list(Intercept = 1, x = 0.5,
                                     sigma = 0.5,
                                     sd_g__Intercept = 0.7))
  d2 <- dd
  d2$y <- s[[1L]]
  f <- frm(form, data = d2)
  expect_lt(abs(sqrt(VarCorr(f)[[1]][1, 1]) - 0.7), 0.25)
  expect_lt(abs(sigma(f) - 0.5), 0.1)
  expect_lt(abs(fixef(f)$mu[["x"]] - 0.5), 0.1)
})

test_that("natural newparams reject what they cannot express", {
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)), y = 0)
  form <- bf(y ~ x + (1 | g)) + gaussian()

  expect_error(frm_simulate(form, dd, newparams = list(zeta = 1)),
               "Unknown newparams")
  # the message lists the vocabulary
  expect_error(frm_simulate(form, dd, newparams = list(zeta = 1)),
               "sd_g__Intercept")
  expect_error(frm_simulate(form, dd, newparams = list(Intercept = 1)),
               "No value for")
  expect_error(frm_simulate(form, dd,
                            newparams = list(Intercept = 1, x = 0.5,
                                             sigma = 1,
                                             sd_g__Intercept = -1)),
               "must be positive")
  expect_error(frm_simulate(form, dd), "needs newparams")

  dd2 <- data.frame(x = rnorm(200), t = factor(rep(1:10, 20)),
                    g = factor(rep(1:20, each = 10)), y = 0)
  expect_error(
    frm_simulate(bf(y ~ x + ar1(t + 0 | g)) + gaussian(), dd2,
                 newparams = list(Intercept = 1)),
    "cannot set the 'ar1' block"
  )
})

test_that("prior draws simulate a prior-predictive sample", {
  set.seed(3)
  dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)), y = 0)
  form <- bf(y ~ x + (1 | g)) + gaussian()
  pl <- set_prior("normal(0, 1)", class = "b") +
    set_prior("normal(0, 2)", class = "Intercept") +
    set_prior("exponential(1)", class = "sd") +
    set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma")

  pp <- frm_simulate(form, dd, priors = pl, nsim = 20, seed = 11)
  pars <- attr(pp, "pars")
  expect_s3_class(pars, "data.frame")
  expect_equal(nrow(pars), 20L)
  expect_true(all(c("Intercept", "x", "sd_g__Intercept",
                    "sigma_Intercept") %in% names(pars)))
  # wide priors: the drawn parameters actually vary
  expect_gt(stats::sd(pars$x), 0.5)
  expect_gt(stats::sd(pars$sd_g__Intercept), 0.2)
  # class "sd" draws live on the natural sd scale
  expect_true(all(pars$sd_g__Intercept > 0))

  # tight priors reproduce the fixed-parameter simulation statistically
  tight <- set_prior("normal(1, 1e-6)", class = "Intercept") +
    set_prior("normal(0.5, 1e-6)", class = "b") +
    set_prior("normal(0.7, 1e-9)", class = "sd") +
    set_prior(paste0("normal(", log(0.6), ", 1e-9)"), class = "Intercept",
              dpar = "sigma")
  a <- frm_simulate(form, dd, priors = tight, nsim = 8, seed = 42)
  b <- frm_simulate(form, dd, nsim = 8, seed = 42,
                    newparams = list(Intercept = 1, x = 0.5, sigma = 0.6,
                                     sd_g__Intercept = 0.7))
  expect_equal(attr(a, "pars")$sd_g__Intercept, rep(0.7, 8),
               tolerance = 1e-6)
  expect_lt(abs(mean(as.matrix(a)) - mean(as.matrix(b))), 0.25)
  expect_lt(abs(stats::sd(as.matrix(a)) - stats::sd(as.matrix(b))), 0.15)
})

test_that("prior draws respect bounds and newparams fill the rest", {
  dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)), y = 0)
  form <- bf(y ~ x + (1 | g)) + gaussian()
  pl <- set_prior("normal(0, 3)", class = "b", lb = 1, ub = 2)
  pp <- frm_simulate(form, dd, priors = pl, nsim = 15, seed = 2,
                     newparams = list(Intercept = 0, x = 0, sigma = 1,
                                      sd_g__Intercept = 0.5))
  drawn <- attr(pp, "pars")$x
  expect_length(drawn, 15L)
  expect_true(all(drawn >= 1 & drawn <= 2))

  # an unreachable bound reports itself rather than spinning
  expect_error(
    frm_simulate(form, dd, nsim = 1, seed = 1,
                 priors = set_prior("normal(0, 0.01)", class = "b",
                                    lb = 100),
                 newparams = list(Intercept = 0, x = 0, sigma = 1,
                                  sd_g__Intercept = 0.5)),
    "did not produce a draw"
  )
})

test_that("simulate() draws latent responses unless censored = TRUE", {
  set.seed(8)
  n <- 200
  dc <- data.frame(x = rnorm(n))
  dc$ylat <- 1 + 0.8 * dc$x + rnorm(n, 0, 1)
  dc$cc <- ifelse(dc$ylat > 1.5, "right", "none")
  dc$y <- pmin(dc$ylat, 1.5)
  fc <- frm(bf(y | cens(cc) ~ x) + gaussian(), data = dc)

  s_lat <- simulate(fc, nsim = 3, seed = 2)
  s_cen <- simulate(fc, nsim = 3, seed = 2, censored = TRUE)
  # the default is the latent response (brms posterior_predict semantics)
  expect_gt(max(s_lat[[1L]]), 1.5)
  # censored = TRUE records every draw inside the observation window
  expect_lte(max(as.matrix(s_cen)), 1.5)
  expect_equal(pmin(as.matrix(s_lat), 1.5), as.matrix(s_cen),
               ignore_attr = TRUE)
  # a censored draw has the same point mass at the censoring point the
  # data have
  expect_gt(sum(s_cen[[1L]] == 1.5), 10L)

  # row-varying censoring times: refused, not silently half-applied
  dv <- dc
  dv$cut <- 1.5 + (seq_len(n) %% 2) * 0.5
  dv$cc <- ifelse(dv$ylat > dv$cut, "right", "none")
  dv$y <- pmin(dv$ylat, dv$cut)
  fv <- frm(bf(y | cens(cc) ~ x) + gaussian(), data = dv)
  expect_error(simulate(fv, nsim = 1, censored = TRUE),
               "type-I censoring")

  # no cens(): the argument has nothing to apply
  fu <- frm(bf(ylat ~ x) + gaussian(), data = dc)
  expect_error(simulate(fu, nsim = 1, censored = TRUE),
               "needs a cens\\(\\) response")
})

test_that("CE prediction intervals respect trials() and trunc()", {
  set.seed(2)
  dn <- data.frame(x = rnorm(150), nt = sample(5:20, 150, TRUE))
  dn$y <- rbinom(150, dn$nt, stats::plogis(0.3 + 0.5 * dn$x))
  fb <- frm(bf(y | trials(nt) ~ x) + binomial(), data = dn)

  # a reference number of trials is meaningless: ask for a real one
  expect_error(
    conditional_effects(fb, effects = "x", method = "predict",
                        resolution = 4, ndraws = 20),
    "trials\\(nt\\)"
  )
  ce <- conditional_effects(fb, effects = "x", method = "predict",
                            resolution = 4, ndraws = 400,
                            conditions = list(nt = 10))
  expect_true(all(ce$x$upper__ <= 10))
  expect_true(all(ce$x$lower__ >= 0))
  # the point estimate is on the count scale the bands live on
  expect_true(all(ce$x$estimate__ > 1))
  expect_true(all(ce$x$estimate__ < 10))

  # literal truncation bounds apply without being pinned
  set.seed(4)
  dt <- data.frame(x = rnorm(300))
  dt$y <- 1 + 0.5 * dt$x + rnorm(300)
  dt <- dt[dt$y > 0.5, ]
  ft <- frm(bf(y | trunc(lb = 0.5) ~ x) + gaussian(), data = dt)
  cet <- conditional_effects(ft, effects = "x", method = "predict",
                             resolution = 4, ndraws = 500)
  expect_true(all(cet$x$lower__ >= 0.5))

  # a variable bound is a real value, so it has to be pinned
  dt$lo <- 0.5
  ft2 <- frm(bf(y | trunc(lb = lo) ~ x) + gaussian(), data = dt)
  expect_error(
    conditional_effects(ft2, effects = "x", method = "predict",
                        resolution = 3, ndraws = 20),
    "trunc\\(lb = lo\\)"
  )
  ce2 <- conditional_effects(ft2, effects = "x", method = "predict",
                             resolution = 3, ndraws = 300,
                             conditions = list(lo = 0.5))
  expect_true(all(ce2$x$lower__ >= 0.5))
})
