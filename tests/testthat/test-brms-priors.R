# Prior placement against brms.
#
# The flat-prior tier next door (test-brms-likelihood.R) proves that
# frmtmb's objective is the same function of the parameters as brms's
# Stan program when neither side carries a prior. This tier carries the
# priors brms's own get_prior() supplies and asks where each one lands.
#
# Two checks, as next door, but the question has changed. With flat
# priors only `adjust_transform = FALSE` could ever match, because
# frmtmb maximizes a likelihood and Stan's Jacobians belong to a
# posterior. With a prior on a TRANSFORMED parameter the Jacobian is
# part of the prior, so:
#
#   A  log_prob under BOTH settings against frmtmb's PENALIZED
#      objective. The difference is reported and attributed, not
#      assumed to be zero.
#   B  grad_log_prob under both settings. The setting whose gradient
#      vanishes names the density frmtmb actually maximized, and that
#      is the answer this file exists to pin.
#
# What it found, and what each test below pins:
#
#   * brms's get_prior() defaults reach frm(prior =) as NOTHING. Every
#     row carries source == "default" and as_priorlist() drops exactly
#     those, so the fit is unpenalized and fit$prior is NULL.
#   * class "sd" is ALREADY brms's placement: the density at the
#     natural sd plus the log-Jacobian. It differs from brms by
#     log(2) per parameter, the half-t renormalizer, which is a
#     constant and moves no mode.
#   * a distributional parameter WITHOUT a linear predictor is the one
#     real placement difference, and frmtmb refuses the row rather than
#     mistranslating it. Its nearest spelling is on the link scale;
#     the natural placement exists internally and reproduces brms
#     exactly.
#   * class "Intercept" is a third thing entirely: same density, same
#     scale, different ARGUMENT, because brms centers its design
#     matrix and frmtmb does not.
#
# Stan compiles here, so the whole file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
# Each shape needs three programs (flat, the honored rows, all the
# defaults), cached under FRMTMB_STAN_CACHE by the same key
# helper-brms.R uses. Four of the six flat programs are byte-identical
# to ones the flat-prior tier already compiles, so a shared cache pays
# for them once.
#
# See dev/brms-priors-findings.md for the tables and the
# recommendation these numbers support.

# ---------------------------------------------------------------------
# The translation surface. No Stan.
# ---------------------------------------------------------------------

test_that("a brms get_prior() table applies no prior at all", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # This is the first thing a user porting a brms script meets, and it
  # is silent apart from one message: the fit succeeds and is
  # unpenalized. Every row get_prior() writes is brms's own default,
  # and as_priorlist() (R/priors.R) drops rows whose source is
  # "default" on purpose, so that frmtmb's defaults and brms's cannot
  # both apply to one parameter. On the frm() path frmtmb has no
  # defaults, so what is left is nothing.
  gp <- brms::get_prior(brms::bf(Reaction ~ Days + (Days | Subject)),
                        data = sleepstudy, family = gaussian())
  expect_true(all(gp$source[nzchar(gp$prior)] == "default"))

  expect_message(
    fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
               data = sleepstudy, prior = gp),
    "brms had filled in as its own defaults")
  expect_null(fit$prior)
  expect_equal(as.numeric(logLik(fit)),
               as.numeric(logLik(frm(bf(Reaction ~ Days +
                                          (Days | Subject)) + gaussian(),
                                     data = sleepstudy))),
               tolerance = 1e-8)
})

test_that("a row the USER edited in a get_prior() table is dropped too", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # The sharpest form of the same defect, and the reason it is a
  # correctness bug rather than an ergonomics gap. `source` records who
  # BUILT the row, not who wrote the density in it, and brms does not
  # update it when a user edits the `prior` cell of a get_prior() table
  # in place, which is the ordinary brms workflow. brms honors such an
  # edit; frmtmb drops it by `source` and then reports it as a row brms
  # filled in itself.
  #
  # Every assertion here pins TODAY's behavior, which is the wrong
  # behavior, so that D1a flips a named expectation rather than
  # arriving unannounced.
  gp <- brms::get_prior(brms::bf(Reaction ~ Days + (1 | Subject)),
                        data = sleepstudy, family = gaussian())
  i <- which(gp$class == "sd" & !nzchar(gp$coef) & !nzchar(gp$group))
  gp$prior[i] <- "normal(0, 20)"
  # brms leaves the row marked as its own default after the edit
  expect_identical(gp$source[[i]], "default")

  m <- bf(Reaction ~ Days + (1 | Subject)) + gaussian()
  # and the message counts the user's own row among brms's defaults
  expect_message(fit <- frm(m, data = sleepstudy, prior = gp),
                 "dropped 3 row")
  expect_null(fit$prior)
  expect_equal(as.numeric(logLik(fit)),
               as.numeric(logLik(frm(m, data = sleepstudy))),
               tolerance = 1e-8)
  # the prior was not vacuous: the same density respelled through
  # set_prior() moves the objective by 1.9 nats, so what the drop costs
  # is the whole of it
  fs <- frm(m, data = sleepstudy,
            prior = set_prior("normal(0, 20)", class = "sd"))
  expect_gt(abs(as.numeric(logLik(fs)) - as.numeric(logLik(fit))), 1)
})

test_that("every default row's fate is one of five, by shape", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # Four distinct refusals, and they are not the same kind of thing.
  # "refused: class"        the class is not one of b/Intercept/sd/cor
  # "refused: distribution" the density is not one of the five parsed
  # "refused: no target"    the class is accepted and addresses nothing
  # "flat slot"             the row is not a prior at all
  #
  # Pinning the whole vector rather than a count means a new brms
  # default, or a new frmtmb class, fails here and names itself.
  st <- function(bform, family, data, frm_model) {
    g <- bp_classify_rows(
      brms::get_prior(bform, data = data, family = family),
      frm(frm_model, data = data))
    g$status[nzchar(g$prior)]
  }

  expect_identical(
    st(brms::bf(Reaction ~ Days + (Days | Subject)), gaussian(),
       sleepstudy, bf(Reaction ~ Days + (Days | Subject)) + gaussian()),
    c("honored", "honored", "honored", "refused: class"))

  set.seed(5)
  do <- data.frame(x = rnorm(300))
  do$y <- ordered(cut(0.9 * do$x + rlogis(300),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))
  # the ordinal threshold prior passes the class gate and then finds
  # nothing: frmtmb's get_prior() offers no slot for a threshold either
  expect_identical(st(brms::bf(y ~ x), brms::cumulative(), do,
                      bf(y ~ x) + cumulative()),
                   "refused: no target")

  set.seed(37)
  dx <- data.frame(x = rnorm(400))
  k <- rbinom(400, 1, 0.35)
  dx$y <- ifelse(k == 1, rnorm(400, 3, 1), rnorm(400, -1, 1))
  # sigma1, sigma2 and brms's mixture theta2 are refused by class; the
  # logistic on Intercept_theta1 is refused because parse_prior_dist()
  # knows five densities and logistic is not one of them
  expect_identical(
    st(brms::bf(y ~ 1, theta1 ~ x),
       brms::mixture(gaussian(), gaussian()), dx,
       bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian())),
    c("refused: class", "refused: class", "refused: class",
      "honored", "honored", "refused: distribution"))

  # The status is right and the ADVICE is not. The special-cased hint
  # tests identical(cls, "theta"), brms spells a mixture proportion
  # theta2, so the generic branch fires and names a spelling that then
  # fails with "Prior target not found". Pinning the text means the
  # one-line fix to that condition has to come past this expectation.
  expect_error(frmtmb:::check_brms_prior_class("theta2", "logistic(0, 1)"),
               'dpar = "theta2"', fixed = TRUE)
})

# ---------------------------------------------------------------------
# The placement identity, on the shape that isolates it.
# ---------------------------------------------------------------------

test_that("row 5: a dpar prior is the whole placement question", {
  skip_unless_brms_fit()

  # The nonlinear shape is where the question has no confounder. Its
  # only default prior is the half-t on sigma; sigma is its only
  # constrained parameter; and brms does not center a nonlinear
  # predictor, so the intercept-centering difference is absent too.
  set.seed(7)
  n <- 120
  dn <- data.frame(x = runif(n, 0, 3))
  dn$y <- 2.5 * exp(-0.8 * dn$x) + rnorm(n, 0, 0.15)
  bform <- brms::bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE)
  r <- bp_shape(bform, gaussian(), dn,
                bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) +
                  gaussian())

  # frm(prior =) honors nothing here, so the "honored" program is the
  # flat one and this row restates the flat tier: frmtmb maximizes the
  # density Stan reports with adjust_transform = FALSE.
  expect_identical(r$rows$status[nzchar(r$rows$prior)], "refused: class")
  expect_lt(abs(r$hon$dF), 1e-8)
  expect_lt(r$hon$gF, 1e-3)
  # and the OTHER setting is off by exactly the derivative of the one
  # log transform in the program, which is 1
  expect_equal(r$hon$gT, 1, tolerance = 1e-3)

  # the hyperparameters come off the row rather than being written in,
  # because brms scales its default by the spread of the response and a
  # hard-coded 2.5 would pin this data set instead of the rule
  h <- bp_hyper(r$rows$prior[[which(nzchar(r$rows$prior))]])

  # (a) the nearest spelling frmtmb's own refusal message suggests,
  #     class = "Intercept" with dpar = "sigma", puts the density on
  #     LOG sigma with no Jacobian
  sl <- as.numeric(r$link$pars[["sigma"]])
  expect_equal(r$link$frm_prior, bp_st(log(sl), h[[1]], h[[2]], h[[3]]),
               tolerance = 1e-10)
  # neither Stan density is the one frmtmb maximized. The two gradients
  # differ by the log transform's derivative, which is 1, so the link
  # fit sits strictly between the two Stan optima rather than at either.
  expect_gt(r$link$gF, 1e-2)
  expect_gt(r$link$gT, 1e-2)

  # (b) the natural placement is brms's, exactly, up to the half-t
  #     renormalizer brms writes and frmtmb does not
  sn <- as.numeric(r$nat$pars[["sigma"]])
  expect_equal(r$nat$frm_prior,
               bp_half_st(sn, h[[1]], h[[2]], h[[3]]) - log(2) + log(sn),
               tolerance = 1e-10)
  # CHECK A: the whole residual is that constant, and the count comes
  # off the program's own lccdf lines rather than being assumed
  expect_identical(bp_half_t_count(r$code$full, r$sdat), 1L)
  expect_equal(r$nat$dT, bp_half_t_const(1), tolerance = 1e-8)
  # CHECK B: and so the AT=TRUE gradient, and only that one, vanishes
  expect_lt(r$nat$gT, 1e-3)
  expect_equal(r$nat$gF, 1, tolerance = 1e-3)
})

test_that("row 5: the natural placement reproduces brms's mode", {
  skip_unless_brms_fit()

  # What a user porting a brms script experiences. Stan's mode under
  # adjust_transform = TRUE is what brms's posterior is a mode of; the
  # two frmtmb spellings are the two things frm(prior =) can be made to
  # say.
  set.seed(7)
  n <- 120
  dn <- data.frame(x = runif(n, 0, 3))
  dn$y <- 2.5 * exp(-0.8 * dn$x) + rnorm(n, 0, 0.15)
  bform <- brms::bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE)
  frm_model <- bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) + gaussian()
  r <- bp_shape(bform, gaussian(), dn, frm_model)

  u <- rstan::unconstrain_pars(r$sf$full, r$nat$pars)
  mode_T <- stats::optim(
    u, function(z) -rstan::log_prob(r$sf$full, z, adjust_transform = TRUE),
    function(z) -rstan::grad_log_prob(r$sf$full, z,
                                      adjust_transform = TRUE),
    method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
  brms_sigma <- as.numeric(
    rstan::constrain_pars(r$sf$full, mode_T$par)[["sigma"]])

  # the natural spelling IS brms's mode
  expect_equal(as.numeric(r$nat$pars[["sigma"]]), brms_sigma,
               tolerance = 1e-6)
  # the link spelling is not, and it lands between the unpenalized
  # estimate and brms's mode rather than to one side of them
  link_sigma <- as.numeric(r$link$pars[["sigma"]])
  mle_sigma <- as.numeric(r$hon$pars[["sigma"]])
  expect_gt(abs(link_sigma - brms_sigma), 1e-5)
  expect_true(link_sigma > mle_sigma && link_sigma < brms_sigma)
})

# ---------------------------------------------------------------------
# The classes that do translate.
# ---------------------------------------------------------------------

test_that("row C: class sd is brms's placement, up to log(2)", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # brms writes a lower-bounded prior as the density MINUS one lccdf
  # per element; frmtmb writes the untruncated density PLUS the log
  # Jacobian. The two therefore differ by log(sd) - log(2) per
  # parameter, and only the second half of that is a constant: the
  # first half is exactly the Jacobian Stan adds under
  # adjust_transform = TRUE. Placing them side by side is what shows
  # that frmtmb is already on brms's side of this question.
  bform <- brms::bf(Reaction ~ Days + (1 | Subject))
  r <- bp_shape(bform, gaussian(), sleepstudy,
                bf(Reaction ~ Days + (1 | Subject)) + gaussian(),
                joint = TRUE)

  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  expect_identical(ent$scale[ent$comp == "theta"], "sd")

  sd1 <- as.numeric(r$hon$pars[["sd_1"]])
  i_sd <- which(r$rows$class == "sd" & nzchar(r$rows$prior))
  h <- bp_hyper(r$rows$prior[[i_sd]])
  expect_equal(ent$value[ent$comp == "theta"],
               bp_half_st(sd1, h[[1]], h[[2]], h[[3]]) - log(2) + log(sd1),
               tolerance = 1e-9)

  # the Intercept row is a DIFFERENT kind of difference: same density,
  # same scale, different argument. brms centers X inside the Stan
  # program, so its Intercept is the intercept at the mean of the
  # predictors and frmtmb's is the intercept at zero.
  i_ic <- which(r$rows$class == "Intercept" & nzchar(r$rows$prior))
  hi <- bp_hyper(r$rows$prior[[i_ic]])
  raw <- fixef(r$fit$hon)$mu[["(Intercept)"]]
  centered <- as.numeric(r$hon$pars[["Intercept"]])
  expect_gt(abs(raw - centered), 1)
  expect_equal(ent$value[ent$comp == "beta"],
               bp_st(raw, hi[[1]], hi[[2]], hi[[3]]), tolerance = 1e-10)

  # CHECK A, decomposed: the residual is the sd Jacobian, minus the
  # half-t renormalizer, plus the centering
  expect_equal(
    r$hon$frm_prior - r$hon$stan_prior_F,
    (log(sd1) - log(2)) +
      (bp_st(raw, hi[[1]], hi[[2]], hi[[3]]) -
         bp_st(centered, hi[[1]], hi[[2]], hi[[3]])),
    tolerance = 1e-8)

  # CHECK B on the z block, which is the only block the joint gradient
  # says anything about (see check C of the flat-prior tier). A prior
  # on the outer parameters leaves the inner problem alone, so the
  # conditional modes are still exactly Stan's.
  expect_lt(r$hon$gFz, 1e-8)
  expect_lt(r$hon$gTz, 1e-8)
  expect_lt(r$nat$gFz, 1e-8)
})

test_that("row C: class cor is the same LKJ in another coordinate", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # brms declares a cholesky_factor_corr and puts lkj_corr_cholesky on
  # it; frmtmb holds a row-normalized Cholesky parameter `t` and
  # carries the SAME density on the correlation matrix onto `t` with
  # that map's exact Jacobian. Both are proper densities on the
  # correlation; they differ by the Jacobian between the two
  # unconstrained coordinates, which is a function of rho and not a
  # constant.
  bform <- brms::bf(Reaction ~ Days + (Days | Subject))
  r <- bp_shape(bform, gaussian(), sleepstudy,
                bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
                joint = TRUE)

  lmat <- r$hon$pars[["L_1"]]
  rho <- lmat[2, 1]
  d <- nrow(lmat)
  eta <- bp_hyper(r$rows$prior[[which(r$rows$class == "cor" &
                                        nzchar(r$rows$prior))]])[[1]]

  # for eta = 1 and d = 2 the LKJ is uniform on rho, so Stan's
  # statement is a constant: log(1/2)
  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  frm_cor <- ent$value[ent$kind == "lkj"]
  # frmtmb's value on its own coordinate, from the closed form of the
  # map: p(t) = p(rho) * |drho/dt| with drho/dt = (1 - rho^2)^(3/2)
  expect_equal(frm_cor,
               -log(2) + (eta + (d - 1) / 2) * log(1 - rho^2),
               tolerance = 1e-9)

  # and the whole residual of check A is still the sum of the named
  # pieces, with the correlation contributing its coordinate change
  sd1 <- as.numeric(r$hon$pars[["sd_1"]])
  i_sd <- which(r$rows$class == "sd" & nzchar(r$rows$prior))
  hs <- bp_hyper(r$rows$prior[[i_sd]])
  i_ic <- which(r$rows$class == "Intercept" & nzchar(r$rows$prior))
  hi <- bp_hyper(r$rows$prior[[i_ic]])
  raw <- fixef(r$fit$hon)$mu[["(Intercept)"]]
  centered <- as.numeric(r$hon$pars[["Intercept"]])
  expect_equal(
    r$hon$frm_prior - r$hon$stan_prior_F,
    sum(log(sd1)) - bp_half_t_const(length(sd1)) +
      (bp_st(raw, hi[[1]], hi[[2]], hi[[3]]) -
         bp_st(centered, hi[[1]], hi[[2]], hi[[3]])) +
      (eta + (d - 1) / 2) * log(1 - rho^2),
    tolerance = 1e-8)
  # the honored program's only lower-bounded prior is the sd vector, so
  # the renormalizer brms writes is one log(2) per standard deviation.
  # Reading the count off the program is what would catch brms changing
  # how it writes a truncated prior.
  expect_identical(bp_half_t_count(r$code$hon, r$sdat),
                   as.integer(length(sd1)))
  expect_equal(sum(bp_half_st(sd1, hs[[1]], hs[[2]], hs[[3]])),
               bp_st(sd1, hs[[1]], hs[[2]], hs[[3]]) +
                 bp_half_t_const(length(sd1)), tolerance = 1e-12)
  expect_lt(r$hon$gFz, 1e-8)
})

test_that("row 1: a dpar WITH a linear predictor already agrees", {
  skip_unless_brms_fit()

  # When sigma has a linear predictor both packages put the prior on
  # the intercept of LOG sigma, so there is no placement question at
  # all and the Jacobian sum over the whole program is zero. What is
  # left is the centering, and with mean-zero predictors that is
  # numerically nothing.
  set.seed(11)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x))
  r <- bp_shape(brms::bf(y ~ x + z, sigma ~ x), gaussian(), dd,
                bf(y ~ x + z, sigma ~ x) + gaussian())

  expect_identical(r$rows$status[nzchar(r$rows$prior)],
                   c("honored", "honored"))
  expect_equal(r$hon$jac, 0, tolerance = 1e-12)
  expect_equal(r$hon$dF, r$hon$dT, tolerance = 1e-12)

  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  expect_identical(sort(unique(ent$scale)), "internal")
  fe <- fixef(r$fit$hon)
  hm <- bp_hyper(r$rows$prior[[which(r$rows$class == "Intercept" &
                                       !nzchar(r$rows$dpar) &
                                       nzchar(r$rows$prior))]])
  hs <- bp_hyper(r$rows$prior[[which(r$rows$class == "Intercept" &
                                       r$rows$dpar == "sigma" &
                                       nzchar(r$rows$prior))]])
  expect_equal(
    r$hon$frm_prior - r$hon$stan_prior_F,
    (bp_st(fe$mu[["(Intercept)"]], hm[[1]], hm[[2]], hm[[3]]) -
       bp_st(as.numeric(r$hon$pars[["Intercept"]]),
             hm[[1]], hm[[2]], hm[[3]])) +
      (bp_st(fe$sigma[["(Intercept)"]], hs[[1]], hs[[2]], hs[[3]]) -
         bp_st(as.numeric(r$hon$pars[["Intercept_sigma"]]),
               hs[[1]], hs[[2]], hs[[3]])),
    tolerance = 1e-10)
})

test_that("S7: the centering is what biases a regression slope", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # The largest effect this tier measures and the only one that moves a
  # regression coefficient rather than a dispersion parameter.
  # Reaction ~ Days with the random effect dropped is the one shape
  # where the centering is structural (mean(Days) = 4.5) AND Stan's
  # optimum is a real mode, so the consequence can be measured instead
  # of argued.
  r <- bp_shape(brms::bf(Reaction ~ Days), gaussian(), sleepstudy,
                bf(Reaction ~ Days) + gaussian())
  expect_identical(r$rows$status[nzchar(r$rows$prior)],
                   c("honored", "refused: class"))

  # The honored program carries brms's Intercept row and nothing else,
  # which is exactly the row frm(prior = ) accepted, so both sides hold
  # the same prior and the only difference left is the argument the
  # density reads.
  mode_T <- function(sf, start) {
    o <- stats::optim(
      start, function(z) -rstan::log_prob(sf, z, adjust_transform = TRUE),
      function(z) -rstan::grad_log_prob(sf, z, adjust_transform = TRUE),
      method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
    rstan::constrain_pars(sf, o$par)
  }
  u <- rstan::unconstrain_pars(r$sf$hon, r$hon$pars)
  brms_days <- mode_T(r$sf$hon, u)[["b"]]
  flat_days <- mode_T(r$sf$flat, u)[["b"]]
  mle_days <- fixef(r$fit0)$mu[["Days"]]
  frm_days <- fixef(r$fit$hon)$mu[["Days"]]
  se_days <- summary(r$fit0)$coefficients$mu["Days", 2]

  # brms's prior does not move the slope: the intercept it constrains
  # is the one at the mean of Days, orthogonal to the slope by
  # construction
  expect_lt(abs(brms_days - flat_days), 1e-3)
  # frmtmb's does, because the intercept IT constrains is the one at
  # Days = 0, which in this design is strongly correlated with the
  # slope. 0.068 standard errors, on 180 rows, from a prior the user
  # believes they carried over unchanged.
  expect_lt(frm_days, mle_days - 0.05)
  expect_equal(abs(frm_days - mle_days) / se_days, 0.0684,
               tolerance = 0.02)

  # and the density difference behind it is the centering term with
  # nothing else on this shape: the two intercepts differ by exactly
  # mean(Days) times the slope
  i_ic <- which(r$rows$class == "Intercept" & nzchar(r$rows$prior))
  hi <- bp_hyper(r$rows$prior[[i_ic]])
  raw <- fixef(r$fit$hon)$mu[["(Intercept)"]]
  centered <- as.numeric(r$hon$pars[["Intercept"]])
  expect_equal(raw + mean(sleepstudy$Days) * frm_days, centered,
               tolerance = 1e-6)
  expect_equal(r$hon$frm_prior - r$hon$stan_prior_F,
               bp_st(raw, hi[[1]], hi[[2]], hi[[3]]) -
                 bp_st(centered, hi[[1]], hi[[2]], hi[[3]]),
               tolerance = 1e-10)
})

test_that("row 12: brms's ordinal threshold prior has no spelling", {
  skip_unless_brms_fit()

  # The one default an ordinal model gets is a student_t on the
  # thresholds. Its class is "Intercept", which the translator
  # accepts, and then resolve_priorlist() finds nothing to attach it
  # to. frmtmb's own get_prior() offers no threshold slot either, so
  # this is a gap in the class vocabulary and not a placement choice.
  set.seed(5)
  n <- 300
  do <- data.frame(x = rnorm(n))
  do$y <- ordered(cut(0.9 * do$x + rlogis(n),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))
  r <- bp_shape(brms::bf(y ~ x), brms::cumulative(), do,
                bf(y ~ x) + cumulative())

  expect_identical(r$rows$status[nzchar(r$rows$prior)],
                   "refused: no target")
  expect_error(
    frm(bf(y ~ x) + cumulative(), data = do,
        prior = bp_frm_prior(r$rows, which(nzchar(r$rows$prior)))),
    "Prior target not found")

  # frmtmb therefore carries no prior at all, and its objective is
  # the flat-prior one: the whole of brms's default is unrepresentable
  expect_equal(r$hon$frm_prior, 0)
  expect_lt(abs(r$hon$dF), 1e-8)
  expect_lt(r$hon$gF, 1e-3)
  # the Jacobian Stan adds for the ordered threshold vector is what
  # separates the two settings, and it is not zero
  expect_gt(abs(r$hon$jac), 1e-3)
})

test_that("row 17: a mixture keeps its intercepts and loses the rest", {
  skip_unless_brms_fit()

  # Four of the mixture's six live defaults are refused. sigma1 and
  # sigma2 are the placement question again, once per component. The
  # logistic on brms's theta is refused for its DENSITY, which no
  # placement decision would fix, and it is refused twice: once on the
  # theta2 class and once on Intercept/theta1.
  set.seed(37)
  n <- 400
  dx <- data.frame(x = rnorm(n))
  k <- rbinom(n, 1, 0.35)
  dx$y <- ifelse(k == 1, rnorm(n, 3, 1), rnorm(n, -1, 1))
  r <- bp_shape(brms::bf(y ~ 1, theta1 ~ x),
                brms::mixture(gaussian(), gaussian()), dx,
                bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()))

  # the two component intercepts translate; brms declares them as one
  # ordered vector, so the prior lands on its entries
  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  expect_identical(ent$comp, c("beta", "beta"))
  expect_identical(sort(unique(ent$scale)), "internal")

  # sigma1 and sigma2 behave exactly as row 5's single sigma does: the
  # natural spelling is brms's up to log(2) PER COMPONENT
  s1 <- as.numeric(r$nat$pars[["sigma1"]])
  s2 <- as.numeric(r$nat$pars[["sigma2"]])
  h <- bp_hyper(r$rows$prior[[which(r$rows$class == "sigma1")]])
  entn <- bp_prior_entries(r$fit$nat, r$prior$nat)
  expect_equal(
    sum(entn$value[entn$comp == "betad"]),
    bp_half_st(c(s1, s2), h[[1]], h[[2]], h[[3]]) -
      bp_half_t_const(2) + log(s1) + log(s2),
    tolerance = 1e-9)
})
