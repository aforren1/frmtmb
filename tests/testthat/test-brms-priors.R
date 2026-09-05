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
# What it now pins, after the PRIORS lane (see dev/priors-findings.md):
#
#   * a brms get_prior() table applies what its rows SAY. Rows used to
#     be dropped by the `source` column, which cost a user their own
#     edited prior; the column is not read any more.
#   * a distributional parameter's own class (`sigma`, `phi`, ...) is
#     ROUTED to the slot holding that parameter's intercept and marked
#     `natural`, so the density is about the parameter, which is where
#     brms puts it. frmtmb's own class = "Intercept" + dpar = spelling
#     still means the LINK scale, and that divergence is measured here
#     too rather than left to a document.
#   * class "Intercept" is evaluated at the intercept brms constrains,
#     the one at the predictor MEANS. The residual that used to be
#     "centering" is zero everywhere in this file.
#   * class "sd" was already brms's placement and still is: the density
#     at the natural sd plus the log-Jacobian, differing from brms by
#     log(2) per lower-bounded element, which is a constant.
#   * an ordinal family's thresholds ARE its class "Intercept", on both
#     sides, with the same map and the same centering.
#
# The consequence, and the reason the checks below are identities: on
# every shape here whose optimum is a real mode, frmtmb's penalized
# objective is brms's posterior density up to one log(2) per
# lower-bounded parameter, so the AT=TRUE gradient vanishes.
#
# Stan compiles here, so the whole file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
# Each shape needs three programs (flat, the honored rows, all the
# defaults), cached under FRMTMB_STAN_CACHE by the same key
# helper-brms.R uses.

# ---------------------------------------------------------------------
# The translation surface. No Stan.
# ---------------------------------------------------------------------

test_that("a brms get_prior() table applies what its rows say", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # This is the first thing a user porting a brms script meets. Every
  # row get_prior() writes carries source == "default", and
  # as_priorlist() used to drop exactly those, so the fit was
  # unpenalized and fit$prior was NULL. It now reads the `prior` string
  # and nothing else, which is brms's own rule.
  gp <- brms::get_prior(brms::bf(Reaction ~ Days + (Days | Subject)),
                        data = sleepstudy, family = gaussian())
  expect_true(all(gp$source[nzchar(gp$prior)] == "default"))

  m <- bf(Reaction ~ Days + (Days | Subject)) + gaussian()
  expect_silent(fit <- frm(m, data = sleepstudy, prior = gp))
  expect_s3_class(fit$prior, "frmtmb_priorlist")
  # one spec per live row, and every one of them resolves
  expect_identical(length(unclass(fit$prior)), sum(nzchar(gp$prior)))
  # and the fit is penalized: the whole table is worth several nats
  expect_gt(abs(as.numeric(logLik(fit)) -
                  as.numeric(logLik(frm(m, data = sleepstudy)))), 1)
})

test_that("a row the USER edited in a get_prior() table is honored", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # The sharpest form of the same defect, and the reason it was a
  # correctness bug rather than an ergonomics gap. `source` records who
  # BUILT the row, not who wrote the density in it, and brms does not
  # update it when a user edits the `prior` cell of a get_prior() table
  # in place, which is the ordinary brms workflow. brms honors such an
  # edit; frmtmb used to drop it by `source` and then report it as a row
  # brms had filled in itself.
  gp <- brms::get_prior(brms::bf(Reaction ~ Days + (1 | Subject)),
                        data = sleepstudy, family = gaussian())
  i <- which(gp$class == "sd" & !nzchar(gp$coef) & !nzchar(gp$group))
  gp$prior[i] <- "normal(0, 20)"
  # brms still leaves the row marked as its own default after the edit,
  # which is exactly why the column cannot be the key
  expect_identical(gp$source[[i]], "default")

  m <- bf(Reaction ~ Days + (1 | Subject)) + gaussian()
  expect_silent(fit <- frm(m, data = sleepstudy, prior = gp))
  # the edited row is in the applied prior, spelled as the user wrote it
  sds <- Filter(function(s) identical(s$class, "sd"), unclass(fit$prior))
  expect_length(sds, 1L)
  expect_equal(sds[[1L]]$dist, prior_normal(0, 20))

  # and it is worth what the user asked for: the entry is the edited
  # density at the natural sd, plus the log-Jacobian class "sd" carries
  e_tab <- bp_prior_entries(fit, fit$prior)
  sd1 <- exp(as.numeric(fit$estimates[["theta"]])[[1L]])
  expect_equal(e_tab$value[e_tab$comp == "theta"],
               stats::dnorm(sd1, 0, 20, log = TRUE) + log(sd1),
               tolerance = 1e-10)
  # and it is the user's density rather than the one brms had written
  # into that cell. normal(0, 20) is tighter than the student_t(3, 0,
  # 59.3) brms puts there, so the standard deviation it shrinks to is
  # smaller: the edit is doing work, in the direction it asks for.
  gd <- brms::get_prior(brms::bf(Reaction ~ Days + (1 | Subject)),
                        data = sleepstudy, family = gaussian())
  fit_def <- frm(m, data = sleepstudy, prior = gd)
  sd_def <- exp(as.numeric(fit_def$estimates[["theta"]])[[1L]])
  expect_lt(sd1, sd_def)
  expect_gt(sd_def - sd1, 0.1)
})

test_that("every default row's fate is one of five, by shape", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # Four distinct refusals, and they are not the same kind of thing.
  # "refused: class"        the class names a structure frmtmb holds
  #                         somewhere else, or not at all
  # "refused: distribution" the density is not one of the seven parsed
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

  # cor, Intercept, sd, sigma: every one of them lands
  expect_identical(
    st(brms::bf(Reaction ~ Days + (Days | Subject)), gaussian(),
       sleepstudy, bf(Reaction ~ Days + (Days | Subject)) + gaussian()),
    c("honored", "honored", "honored", "honored"))

  set.seed(5)
  do <- data.frame(x = rnorm(300))
  do$y <- ordered(cut(0.9 * do$x + rlogis(300),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))
  # an ordinal family's thresholds ARE its class "Intercept", here as in
  # brms, so the row that used to find no target now lands on them
  expect_identical(st(brms::bf(y ~ x), brms::cumulative(), do,
                      bf(y ~ x) + cumulative()),
                   "honored")

  set.seed(37)
  dx <- data.frame(x = rnorm(400))
  k <- rbinom(400, 1, 0.35)
  dx$y <- ifelse(k == 1, rnorm(400, 3, 1), rnorm(400, -1, 1))
  # sigma1 and sigma2 are routed to their own parameters, and the
  # logistic on Intercept_theta1 parses now. What is left refused is
  # brms's theta2: the word names frmtmb's raw covariance vector, and
  # the component brms holds as its reference is a parameter on neither
  # side. brms writes no Stan statement for it either.
  expect_identical(
    st(brms::bf(y ~ 1, theta1 ~ x),
       brms::mixture(gaussian(), gaussian()), dx,
       bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian())),
    c("honored", "honored", "refused: class",
      "honored", "honored", "honored"))

  # The refusal is right and its ADVICE now is too. It used to name
  # class = "Intercept", dpar = "theta2", which then failed with "Prior
  # target not found", because that component has no linear predictor.
  msg <- tryCatch(frmtmb:::check_brms_prior_class("theta2",
                                                  "logistic(0, 1)"),
                  error = conditionMessage)
  expect_match(msg, "mixture proportion", fixed = TRUE)
  expect_match(msg, "reference component", fixed = TRUE)
  expect_false(grepl('dpar = "theta2"', msg, fixed = TRUE))

  # and the classes that name a structure frmtmb keeps elsewhere say
  # where it is rather than offering a spelling that does not fit
  expect_error(frmtmb:::check_brms_prior_class("sds",
                                               "student_t(3, 0, 1)"),
               "class = \"sd\" with group", fixed = TRUE)
})

# ---------------------------------------------------------------------
# The placement identity, on the shape that isolates it.
# ---------------------------------------------------------------------

test_that("row 5: a dpar prior lands on the parameter brms means", {
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

  # the row is honored, and honored on the NATURAL scale: the entry
  # carries the same change of variables class "sd" uses, because
  # sigma's link is the log
  expect_identical(r$rows$status[nzchar(r$rows$prior)], "honored")
  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  expect_identical(ent$scale, "sd")
  expect_identical(ent$comp, "betad")

  # the hyperparameters come off the row rather than being written in,
  # because brms scales its default by the spread of the response and a
  # hard-coded 2.5 would pin this data set instead of the rule
  h <- bp_hyper(r$rows$prior[[which(nzchar(r$rows$prior))]])

  # (a) frmtmb's density IS brms's, up to the half-t renormalizer brms
  #     writes and frmtmb does not
  sn <- as.numeric(r$hon$pars[["sigma"]])
  expect_equal(r$hon$frm_prior,
               bp_half_st(sn, h[[1]], h[[2]], h[[3]]) - log(2) + log(sn),
               tolerance = 1e-10)
  # CHECK A: the whole residual is that constant, and the count comes
  # off the program's own lccdf lines rather than being assumed
  expect_identical(bp_half_t_count(r$code$hon, r$sdat), 1L)
  expect_equal(r$hon$dT, bp_half_t_const(1), tolerance = 1e-8)
  # CHECK B: and so the AT=TRUE gradient, and only that one, vanishes.
  # frmtmb maximizes exactly the density brms samples.
  expect_lt(r$hon$gT, 1e-3)
  expect_equal(r$hon$gF, 1, tolerance = 1e-3)

  # (b) frmtmb's OWN spelling for the same parameter still means LOG
  #     sigma. It is left alone deliberately: flipping it would change
  #     what an existing frmtmb script means (?set_prior records the
  #     divergence). Measured here so its size cannot rot.
  sl <- as.numeric(r$link$pars[["sigma"]])
  expect_equal(r$link$frm_prior, bp_st(log(sl), h[[1]], h[[2]], h[[3]]),
               tolerance = 1e-10)
  # neither Stan density is the one that spelling maximized: the two
  # gradients differ by the log transform's derivative, which is 1, so
  # the link fit sits strictly between the two Stan optima
  expect_gt(r$link$gF, 1e-2)
  expect_gt(r$link$gT, 1e-2)
})

test_that("row 5: the translated table reproduces brms's mode", {
  skip_unless_brms_fit()

  # What a user porting a brms script experiences. Stan's mode under
  # adjust_transform = TRUE is what brms's posterior is a mode of.
  set.seed(7)
  n <- 120
  dn <- data.frame(x = runif(n, 0, 3))
  dn$y <- 2.5 * exp(-0.8 * dn$x) + rnorm(n, 0, 0.15)
  bform <- brms::bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE)
  frm_model <- bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) + gaussian()
  r <- bp_shape(bform, gaussian(), dn, frm_model)

  u <- rstan::unconstrain_pars(r$sf$full, r$hon$pars)
  mode_T <- stats::optim(
    u, function(z) -rstan::log_prob(r$sf$full, z, adjust_transform = TRUE),
    function(z) -rstan::grad_log_prob(r$sf$full, z,
                                      adjust_transform = TRUE),
    method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
  brms_sigma <- as.numeric(
    rstan::constrain_pars(r$sf$full, mode_T$par)[["sigma"]])

  # the translated table IS brms's mode
  expect_equal(as.numeric(r$hon$pars[["sigma"]]), brms_sigma,
               tolerance = 1e-6)
  # frmtmb's own link spelling is not, and it lands between the
  # unpenalized estimate and brms's mode rather than to one side
  link_sigma <- as.numeric(r$link$pars[["sigma"]])
  mle_sigma <- as.numeric(exp(r$fit0$estimates$betad))
  expect_gt(abs(link_sigma - brms_sigma), 1e-5)
  expect_true(link_sigma > mle_sigma && link_sigma < brms_sigma)
})

# ---------------------------------------------------------------------
# The classes that translate.
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
  # adjust_transform = TRUE.
  bform <- brms::bf(Reaction ~ Days + (1 | Subject))
  r <- bp_shape(bform, gaussian(), sleepstudy,
                bf(Reaction ~ Days + (1 | Subject)) + gaussian(),
                joint = TRUE)

  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  expect_identical(ent$scale[ent$comp == "theta"], "sd")
  # sigma's own row rides on the same change of variables
  expect_identical(ent$scale[ent$comp == "betad"], "sd")

  sd1 <- as.numeric(r$hon$pars[["sd_1"]])
  sig <- as.numeric(r$hon$pars[["sigma"]])
  i_sd <- which(r$rows$class == "sd" & nzchar(r$rows$prior))
  h <- bp_hyper(r$rows$prior[[i_sd]])
  expect_equal(ent$value[ent$comp == "theta"],
               bp_half_st(sd1, h[[1]], h[[2]], h[[3]]) - log(2) + log(sd1),
               tolerance = 1e-9)

  # the Intercept row is on the intercept brms constrains: the one at
  # the MEAN of the predictors, not the one at zero. The entry says so,
  # and the two arguments really are different numbers here, because
  # mean(Days) is 4.5.
  i_ic <- which(r$rows$class == "Intercept" & nzchar(r$rows$prior))
  hi <- bp_hyper(r$rows$prior[[i_ic]])
  raw <- fixef(r$fit$hon)$mu[["(Intercept)"]]
  centered <- as.numeric(r$hon$pars[["Intercept"]])
  expect_gt(abs(raw - centered), 1)
  expect_true(ent$centered[ent$comp == "beta"])
  expect_equal(ent$value[ent$comp == "beta"],
               bp_st(centered, hi[[1]], hi[[2]], hi[[3]]),
               tolerance = 1e-8)

  # CHECK A, decomposed: the residual is the two Jacobians minus one
  # half-t renormalizer each, and NOTHING from the intercept
  expect_identical(bp_half_t_count(r$code$hon, r$sdat), 2L)
  expect_equal(
    r$hon$frm_prior - r$hon$stan_prior_F,
    log(sd1) + log(sig) - bp_half_t_const(2),
    tolerance = 1e-6)
  expect_equal(r$hon$dT, bp_half_t_const(2), tolerance = 1e-6)

  # CHECK B on the z block, which is the only block the joint gradient
  # says anything about (see check C of the flat-prior tier). A prior
  # on the outer parameters leaves the inner problem alone, so the
  # conditional modes are still exactly Stan's.
  expect_lt(r$hon$gFz, 1e-8)
  expect_lt(r$hon$gTz, 1e-8)
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
  # constant. That is the one residual this shape keeps.
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

  # and the whole residual of check A is the sum of the named pieces,
  # with no centering term left in it
  sd1 <- as.numeric(r$hon$pars[["sd_1"]])
  sig <- as.numeric(r$hon$pars[["sigma"]])
  expect_equal(
    r$hon$frm_prior - r$hon$stan_prior_F,
    sum(log(sd1)) + log(sig) - bp_half_t_const(length(sd1) + 1L) +
      (eta + (d - 1) / 2) * log(1 - rho^2),
    tolerance = 1e-6)
  # every lower-bounded prior in the honored program is one of those
  # standard deviations or sigma, so the renormalizer count is theirs.
  # Reading it off the program is what would catch brms changing how it
  # writes a truncated prior.
  expect_identical(bp_half_t_count(r$code$hon, r$sdat),
                   as.integer(length(sd1) + 1L))
  hs <- bp_hyper(r$rows$prior[[which(r$rows$class == "sd" &
                                       nzchar(r$rows$prior))]])
  expect_equal(sum(bp_half_st(sd1, hs[[1]], hs[[2]], hs[[3]])),
               bp_st(sd1, hs[[1]], hs[[2]], hs[[3]]) +
                 bp_half_t_const(length(sd1)), tolerance = 1e-12)
  expect_lt(r$hon$gFz, 1e-8)
})

test_that("row 1: a dpar WITH a linear predictor agrees exactly", {
  skip_unless_brms_fit()

  # When sigma has a linear predictor both packages put the prior on
  # the intercept of LOG sigma, so there is no placement question and
  # the Jacobian sum over the whole program is zero. What used to be
  # left was the centering, and both intercepts are centered now, so
  # the residual is zero rather than merely small.
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
  # both intercepts carry the centering offset, one per sub-formula,
  # which is how brms writes them (means_X and means_X_sigma)
  expect_true(all(ent$centered))

  # nothing at all is left over: frmtmb's log prior IS the Stan
  # program's, to machine precision
  expect_equal(r$hon$frm_prior - r$hon$stan_prior_F, 0,
               tolerance = 1e-10)
  hm <- bp_hyper(r$rows$prior[[which(r$rows$class == "Intercept" &
                                       !nzchar(r$rows$dpar) &
                                       nzchar(r$rows$prior))]])
  hs <- bp_hyper(r$rows$prior[[which(r$rows$class == "Intercept" &
                                       r$rows$dpar == "sigma" &
                                       nzchar(r$rows$prior))]])
  expect_equal(
    r$hon$frm_prior,
    bp_st(as.numeric(r$hon$pars[["Intercept"]]),
          hm[[1]], hm[[2]], hm[[3]]) +
      bp_st(as.numeric(r$hon$pars[["Intercept_sigma"]]),
            hs[[1]], hs[[2]], hs[[3]]),
    tolerance = 1e-10)
})

test_that("S7: the centering no longer biases a regression slope", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  # The largest effect this tier ever measured, and the only one that
  # moved a regression coefficient rather than a dispersion parameter.
  # Reaction ~ Days with the random effect dropped is the one shape
  # where the centering is structural (mean(Days) = 4.5) AND Stan's
  # optimum is a real mode, so the consequence is measured rather than
  # argued. It used to be 0.0684 standard errors of slope bias.
  r <- bp_shape(brms::bf(Reaction ~ Days), gaussian(), sleepstudy,
                bf(Reaction ~ Days) + gaussian())
  expect_identical(r$rows$status[nzchar(r$rows$prior)],
                   c("honored", "honored"))

  mode_T <- function(sf, start) {
    o <- stats::optim(
      start, function(z) -rstan::log_prob(sf, z, adjust_transform = TRUE),
      function(z) -rstan::grad_log_prob(sf, z, adjust_transform = TRUE),
      method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
    rstan::constrain_pars(sf, o$par)
  }
  u <- rstan::unconstrain_pars(r$sf$hon, r$hon$pars)
  brms_mode <- mode_T(r$sf$hon, u)
  mle_days <- fixef(r$fit0)$mu[["Days"]]
  frm_days <- fixef(r$fit$hon)$mu[["Days"]]
  se_days <- summary(r$fit0)$coefficients$mu["Days", 2]

  # brms's prior does not move the slope: the intercept it constrains
  # is the one at the mean of Days, orthogonal to the slope by
  # construction. frmtmb's prior is now on the same intercept, so it
  # does not move the slope either.
  expect_lt(abs(frm_days - mle_days) / se_days, 0.005)
  # and frmtmb lands on brms's mode, in every parameter at once
  expect_equal(frm_days, as.numeric(brms_mode[["b"]]), tolerance = 1e-5)
  expect_equal(as.numeric(r$hon$pars[["Intercept"]]),
               as.numeric(brms_mode[["Intercept"]]), tolerance = 1e-5)
  expect_equal(as.numeric(r$hon$pars[["sigma"]]),
               as.numeric(brms_mode[["sigma"]]), tolerance = 1e-5)

  # the identity behind it: the argument the density reads is the raw
  # intercept plus mean(Days) times the slope
  raw <- fixef(r$fit$hon)$mu[["(Intercept)"]]
  centered <- as.numeric(r$hon$pars[["Intercept"]])
  expect_equal(raw + mean(sleepstudy$Days) * frm_days, centered,
               tolerance = 1e-6)
  # so the whole check-A residual is sigma's Jacobian and its
  # renormalizer, with nothing from the intercept
  sig <- as.numeric(r$hon$pars[["sigma"]])
  expect_equal(r$hon$frm_prior - r$hon$stan_prior_F,
               log(sig) - bp_half_t_const(1), tolerance = 1e-8)
  expect_equal(r$hon$dT, bp_half_t_const(1), tolerance = 1e-8)
  expect_lt(r$hon$gT, 1e-3)
})

test_that("row 12: brms's ordinal threshold prior lands on tau_raw", {
  skip_unless_brms_fit()

  # The one default an ordinal model gets is a student_t on the
  # thresholds, and its class is "Intercept" because in brms the
  # thresholds ARE the intercept. frmtmb holds them in `tau_raw` as
  # (tau_1, log increments), which is the same map Stan's `ordered`
  # type applies, and the design is centered on both sides. So the two
  # densities agree exactly, with no renormalizer to separate them.
  set.seed(5)
  n <- 300
  do <- data.frame(x = rnorm(n))
  do$y <- ordered(cut(0.9 * do$x + rlogis(n),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))
  r <- bp_shape(brms::bf(y ~ x), brms::cumulative(), do,
                bf(y ~ x) + cumulative())

  expect_identical(r$rows$status[nzchar(r$rows$prior)], "honored")
  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  expect_identical(ent$comp, "tau_raw")
  expect_identical(ent$scale, "ordthres")
  expect_true(ent$centered)

  # the thresholds the density reads are Stan's `Intercept` vector,
  # which is the frmtmb threshold vector minus mean(X) times the slope
  raw <- r$fit$hon$estimates[["tau_raw"]]
  tau <- c(raw[1], raw[1] + cumsum(exp(raw[-1])))
  slope <- fixef(r$fit$hon)$mu[["x"]]
  expect_equal(tau - mean(do$x) * slope,
               as.numeric(r$hon$pars[["Intercept"]]), tolerance = 1e-6)

  # and the density itself is brms's, with the ordered map's Jacobian
  h <- bp_hyper(r$rows$prior[[which(nzchar(r$rows$prior))]])
  expect_equal(r$hon$frm_prior,
               bp_st(as.numeric(r$hon$pars[["Intercept"]]),
                     h[[1]], h[[2]], h[[3]]) + sum(raw[-1]),
               tolerance = 1e-9)
  # CHECK A and B: no lower-bounded parameter anywhere, so there is no
  # renormalizer and the residual is zero, not a constant
  expect_identical(bp_half_t_count(r$code$hon, r$sdat), 0L)
  expect_equal(r$hon$dT, 0, tolerance = 1e-8)
  expect_lt(r$hon$gT, 1e-2)
  expect_equal(r$hon$gF, 1, tolerance = 1e-2)

  # the internal route still works and still means the internal scale,
  # which is the escape hatch ?set_prior documents
  # one entry per threshold, and no map applied to any of them
  ri <- resolve_prior_input(r$fit0, list(tau_raw = prior_normal(0, 5)))
  expect_length(ri$entries, length(raw))
  expect_identical(unique(vapply(ri$entries, function(e) e$scale, "")),
                   "internal")
})

test_that("row 17: a mixture keeps everything but its reference theta", {
  skip_unless_brms_fit()

  # Five of the mixture's six live defaults land. sigma1 and sigma2 are
  # the placement question again, once per component, and the logistic
  # on brms's modeled theta parses now. What is left refused is the
  # class "theta2" row, which brms writes no Stan statement for either.
  set.seed(37)
  n <- 400
  dx <- data.frame(x = rnorm(n))
  k <- rbinom(n, 1, 0.35)
  dx$y <- ifelse(k == 1, rnorm(n, 3, 1), rnorm(n, -1, 1))
  r <- bp_shape(brms::bf(y ~ 1, theta1 ~ x),
                brms::mixture(gaussian(), gaussian()), dx,
                bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()))

  ent <- bp_prior_entries(r$fit$hon, r$prior$hon)
  # two component sigmas on their own scale, two component intercepts
  # on theirs, and the modeled proportion's intercept
  expect_identical(sort(ent$scale),
                   c("internal", "internal", "internal", "sd", "sd"))
  expect_true("logistic" %in% ent$kind)
  # y ~ 1 gives the component intercepts no predictor to be centered
  # against; theta1 ~ x gives its intercept one
  expect_identical(ent$centered, c(FALSE, FALSE, FALSE, FALSE, TRUE))

  # the two sigmas behave exactly as row 5's single sigma does: brms's
  # density plus the Jacobian, minus one renormalizer per component
  s1 <- as.numeric(r$hon$pars[["sigma1"]])
  s2 <- as.numeric(r$hon$pars[["sigma2"]])
  h <- bp_hyper(r$rows$prior[[which(r$rows$class == "sigma1")]])
  expect_equal(
    sum(ent$value[ent$scale == "sd"]),
    bp_half_st(c(s1, s2), h[[1]], h[[2]], h[[3]]) -
      bp_half_t_const(2) + log(s1) + log(s2),
    tolerance = 1e-9)

  # and the whole check-A residual is those two Jacobians and nothing
  # else: the intercepts and the logistic agree exactly
  expect_identical(bp_half_t_count(r$code$hon, r$sdat), 2L)
  expect_equal(r$hon$frm_prior - r$hon$stan_prior_F,
               log(s1) + log(s2) - bp_half_t_const(2),
               tolerance = 1e-6)
})
