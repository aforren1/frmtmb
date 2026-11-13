# Draws from frm_sample(laplace = TRUE) hold the outer parameters alone:
# the random effects, the smooth coefficients and the mi() values are
# integrated out. The draws methods read such draws in their own layout,
# refuse by name a quantity that needs what was integrated out, and
# compute the rest exactly. Before, every method read them in the full
# layout: posterior_epred() returned NaN (4500 of 4500 cells, seed 1212
# in dev/arcovsample-rev-11-laplace.R) and, at newdata, finite values
# whose group effects were theta_1 and lp__ (dev/sampfix-06-probes.R).
# The pairs below share every outer value, so a quantity the laplace
# draws can compute must equal the full draws' one bit for bit.

lap_data <- function(seed = 1212L) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(1:6, each = 5L)), t = rep(1:5, 6L))
  dd$x <- stats::rnorm(nrow(dd))
  dd$y <- 0.5 + 0.4 * dd$x + stats::rnorm(6L, 0, 0.5)[as.integer(dd$g)] +
    stats::rnorm(nrow(dd), 0, 0.7)
  dd
}

# full draws around the ML estimates, and the same draws with every
# integrated column removed, which is the layout laplace draws have
lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  stopifnot(length(est) == length(lab))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  outer <- setdiff(colnames(M), inner)
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, outer, drop = FALSE]
  list(full = full, lap = lap)
}

test_that("the predictive methods refuse laplace draws that need b", {
  dd <- lap_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  p <- lap_pair(fit)
  expect_true(frmtmb.sample:::draws_is_laplace(p$lap))
  nd <- data.frame(x = c(-1, 1), g = factor(c(1, 2), levels = 1:6))
  msg <- "frm_sample(laplace = TRUE), which integrates them out"
  expect_error(posterior_epred(p$lap), msg, fixed = TRUE)
  expect_error(posterior_epred(p$lap), "posterior_epred() needs",
               fixed = TRUE)
  expect_error(posterior_linpred(p$lap), msg, fixed = TRUE)
  expect_error(posterior_predict(p$lap), msg, fixed = TRUE)
  # the case that was finite and wrong: a group effect read off theta_1
  expect_error(posterior_epred(p$lap, newdata = nd), msg, fixed = TRUE)
  expect_error(posterior_predict(p$lap, newdata = nd), msg, fixed = TRUE)
  expect_error(fitted(p$lap), msg, fixed = TRUE)
  expect_error(predict(p$lap), msg, fixed = TRUE)
  expect_error(residuals(p$lap), msg, fixed = TRUE)
  expect_error(predictive_interval(p$lap), msg, fixed = TRUE)
  expect_error(bayes_R2(p$lap), msg, fixed = TRUE)
  expect_error(suppressMessages(pp_check(p$lap, ndraws = 2)), msg,
               fixed = TRUE)
  # the refusal points at the population-level quantity it can compute
  expect_error(posterior_epred(p$lap), "re_formula = NA", fixed = TRUE)
})

test_that("what laplace draws can compute equals the full draws' answer", {
  dd <- lap_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  p <- lap_pair(fit)
  nd <- data.frame(x = c(-1, 1), g = factor(c(1, 2), levels = 1:6))
  expect_identical(posterior_epred(p$lap, re_formula = NA),
                   posterior_epred(p$full, re_formula = NA))
  expect_identical(posterior_epred(p$lap, newdata = nd, re_formula = NA),
                   posterior_epred(p$full, newdata = nd, re_formula = NA))
  expect_identical(posterior_linpred(p$lap, re_formula = NA),
                   posterior_linpred(p$full, re_formula = NA))
  # sigma has no group-level term, so it needs no b at any re_formula
  expect_identical(posterior_epred(p$lap, dpar = "sigma"),
                   posterior_epred(p$full, dpar = "sigma"))
  # the predictive draw reads sigma, which the full layout misplaced
  set.seed(4)
  a <- posterior_predict(p$lap, re_formula = NA)
  set.seed(4)
  b <- posterior_predict(p$full, re_formula = NA)
  expect_identical(a, b)
  # a variance component is an outer parameter
  expect_identical(
    hypothesis(p$lap, "sd_g__Intercept > 0", class = NULL)$hypothesis,
    hypothesis(p$full, "sd_g__Intercept > 0", class = NULL)$hypothesis)
  expect_identical(VarCorr(p$lap), VarCorr(p$full))
  # conditional_effects() drops the group effects by default, and at
  # re_formula = NULL draws a new level from each draw's own theta; the
  # observed levels' b is read by neither
  expect_identical(conditional_effects(p$lap, effects = "x",
                                       resolution = 5),
                   conditional_effects(p$full, effects = "x",
                                       resolution = 5))
  expect_identical(conditional_effects(p$lap, effects = "x",
                                       resolution = 5, re_formula = NULL,
                                       seed = 2),
                   conditional_effects(p$full, effects = "x",
                                       resolution = 5, re_formula = NULL,
                                       seed = 2))
})

test_that("the refusal leaves the caller's random-number stream alone", {
  dd <- lap_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  p <- lap_pair(fit)
  set.seed(8)
  posterior_epred(p$lap, re_formula = NA)
  a <- stats::runif(1)
  set.seed(8)
  b <- stats::runif(1)
  expect_identical(a, b)
})

test_that("a smooth keeps the refusal at re_formula = NA", {
  dd <- lap_data()
  fit <- suppressWarnings(frm(bf(y ~ s(x, k = 5)), family = gaussian(),
                              data = dd))
  p <- lap_pair(fit)
  expect_true(frmtmb.sample:::draws_is_laplace(p$lap))
  expect_error(posterior_epred(p$lap, re_formula = NA),
               "integrates them out", fixed = TRUE)
  expect_error(conditional_effects(p$lap, effects = "x", resolution = 5),
               "conditional_effects() needs the random effects", fixed = TRUE)
  expect_identical(posterior_epred(p$lap, dpar = "sigma"),
                   posterior_epred(p$full, dpar = "sigma"))
})

test_that("a curve at re_formula = NULL that reads a smooth is refused", {
  # re_formula = NULL draws a new group level per draw, and the smooth's
  # coefficients are still read: the refusal must see that, not an
  # error of its own that it takes for one unrelated to the fill (R CMD
  # check found such a probe calling a helper that no longer existed,
  # and the curves then came back NA without a word)
  dd <- lap_data()
  fit <- suppressWarnings(frm(bf(y ~ s(x, k = 5) + (1 | g)),
                              family = gaussian(), data = dd))
  p <- lap_pair(fit)
  expect_error(conditional_effects(p$lap, effects = "x", resolution = 5,
                                   re_formula = NULL, seed = 2),
               "conditional_effects() needs the random effects",
               fixed = TRUE)
})

test_that("mi() values are integrated too, and only their reader refuses", {
  dd <- lap_data()
  dd$x[c(3, 11, 19)] <- NA
  fit <- frm(bf(y ~ mi(x)) + bf(x | mi() ~ 1) + set_rescor(FALSE),
             family = gaussian(), data = dd)
  expect_false(length(fit$frame[["re_blocks"]]) > 0L)
  p <- lap_pair(fit)
  expect_error(posterior_epred(p$lap, resp = "y"), "integrates them out",
               fixed = TRUE)
  expect_identical(posterior_epred(p$lap, resp = "x"),
                   posterior_epred(p$full, resp = "x"))
  expect_error(log_lik(p$lap), "integrates them out", fixed = TRUE)
})

test_that("pp_mixture() refuses laplace draws of a grouped mixture", {
  set.seed(6)
  dd <- data.frame(g = factor(rep(1:8, each = 15L)))
  dd$y <- c(stats::rnorm(60, -2), stats::rnorm(60, 3)) +
    stats::rnorm(8, 0, 0.3)[dd$g]
  fit <- suppressWarnings(frm(bf(y ~ 1 + (1 | g)),
                              family = frmtmb::mixture(gaussian(),
                                                       gaussian()),
                              data = dd))
  p <- lap_pair(fit)
  expect_error(pp_mixture(p$lap), "pp_mixture() needs", fixed = TRUE)
  # pp_mixture() takes no re_formula, so the refusal does not offer it
  e <- tryCatch(pp_mixture(p$lap), error = function(e) e)
  expect_no_match(conditionMessage(e), "re_formula", fixed = TRUE)
})

test_that("a later draw that overflows on its own is not taken for a read", {
  # exp(720) is Inf at draw 3 for a reason that has nothing to do with
  # the integrated values, and the full draws return it; a read of those
  # values comes out NA, never Inf (dev/sampfix-12-fillpaths.R)
  set.seed(77)
  dd <- data.frame(g = factor(rep(1:8, each = 10)))
  dd$x <- stats::rnorm(nrow(dd))
  u <- stats::rnorm(8, 0, 0.6)[dd$g]
  dd$pos <- exp(0.2 + 0.1 * dd$x + 0.3 * u + stats::rnorm(80, 0, 0.2))
  fit <- suppressWarnings(frm(bf(pos ~ x + (1 | g)), family = lognormal(),
                              data = dd))
  p <- lap_pair(fit)
  p$full$draws[3L, "b_Intercept"] <- 720
  p$lap$draws[3L, "b_Intercept"] <- 720
  a <- posterior_epred(p$lap, re_formula = NA)
  expect_true(any(is.infinite(a[3L, ])))
  expect_identical(a, posterior_epred(p$full, re_formula = NA))
})

test_that("a refusal names the function the user called", {
  dd <- lap_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  p <- lap_pair(fit)
  starts <- function(expr, what) {
    e <- tryCatch(suppressMessages(expr), error = function(e) e)
    expect_s3_class(e, "error")
    expect_true(startsWith(conditionMessage(e), what),
                info = conditionMessage(e))
  }
  # each used to name the inner function it computes through
  starts(fitted(p$lap), "fitted() needs")
  starts(fitted(p$lap, scale = "linear"), "fitted() needs")
  starts(predict(p$lap), "predict() needs")
  starts(residuals(p$lap), "residuals() needs")
  starts(predictive_error(p$lap), "predictive_error() needs")
  starts(predictive_interval(p$lap), "predictive_interval() needs")
  starts(bayes_R2(p$lap), "bayes_R2() needs")
  starts(pp_check(p$lap, ndraws = 2), "pp_check() needs")
  starts(loo(p$lap), "loo() needs")
  starts(waic(p$lap), "waic() needs")
  starts(psis(p$lap), "psis() needs")
  starts(loo_compare(p$lap, p$lap), "loo_compare() needs")
  starts(hypothesis(p$lap, "Intercept > 0", scope = "ranef", group = "g"),
         "hypothesis() has no draws")
  starts(coef(p$lap), "coef() has no draws")
  starts(ranef(p$lap), "ranef() has no draws")
  starts(log_lik(p$lap), "log_lik() needs")
  # and the name does not stick to a later direct call
  starts(posterior_epred(p$lap), "posterior_epred() needs")
})

test_that("re_formula = NA is suggested only where it would compute", {
  dd <- lap_data()
  hint <- "re_formula = NA leaves the group-level effects out"
  msg <- function(expr) {
    conditionMessage(tryCatch(suppressMessages(expr), error = function(e) e))
  }
  p <- lap_pair(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
  expect_match(msg(posterior_epred(p$lap)), hint, fixed = TRUE)
  expect_match(msg(fitted(p$lap)), hint, fixed = TRUE)
  # bayes_R2() takes brms's re_formula since lane postfix, and with NA
  # it computes from these draws, so the hint is right there too
  expect_match(msg(bayes_R2(p$lap)), hint, fixed = TRUE)
  expect_match(msg(bayes_R2(p$lap)), "without laplace = TRUE",
               fixed = TRUE)
  expect_true(is.matrix(suppressMessages(bayes_R2(p$lap, re_formula = NA))))
  # a smooth is read at re_formula = NA too, so NA would be refused again
  ps <- lap_pair(suppressWarnings(frm(bf(y ~ s(x, k = 5) + (1 | g)),
                                      family = gaussian(), data = dd)))
  expect_no_match(msg(posterior_epred(ps$lap)), "re_formula", fixed = TRUE)
  expect_error(posterior_epred(ps$lap, re_formula = NA),
               "integrates them out", fixed = TRUE)
})

test_that("a draw that reads an integrated value after the first is refused", {
  # The probe runs at the first draw. With k = 0 there, exp(a)^k is
  # NA^0 = 1 and the probe sees no read, while every later draw read `a`'s
  # group effect and was returned as NaN (reviewer's construction,
  # dev/sampfix-rev-01-probe.R, data seed 77)
  set.seed(77)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- stats::rnorm(nrow(dn))
  u <- stats::rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 +
    stats::rnorm(nrow(dn), 0, 0.3)
  fnl <- suppressWarnings(suppressMessages(
    frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
           nl = TRUE), family = gaussian(), data = dn)))
  p <- lap_pair(fnl)
  p$lap$draws[1L, "b_k_Intercept"] <- 0
  expect_error(posterior_epred(p$lap), "posterior_epred() needs",
               fixed = TRUE)
  expect_error(posterior_predict(p$lap), "posterior_predict() needs",
               fixed = TRUE)
  # the case the probe catches on its own stays caught
  expect_error(posterior_epred(p$lap, draw_ids = 2:6),
               "integrates them out", fixed = TRUE)
  # and a result the integrated values do not reach is still returned
  expect_identical(posterior_epred(p$lap, dpar = "sigma"),
                   {
                     full <- p$full
                     full$draws[1L, "b_k_Intercept"] <- 0
                     posterior_epred(full, dpar = "sigma")
                   })
})

test_that("a first draw NaN for its own reason does not hide a later read", {
  # Backlog 1456 (dev/sampfix-rev-r2-01-watch.R, data seed 77): with
  # c1 = -1 and k = 0 on draw 1, log(c1) is NaN at every fill, so the
  # watch took draw 1's all-NA pattern as its reference and returned 480
  # of 480 cells NaN where the full draws have 80. Each draw is now
  # compared with itself at two finite fills.
  set.seed(77)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- stats::rnorm(nrow(dn))
  u <- stats::rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 0.5 * dn$x + exp(0.3 + u)^0.8 + stats::rnorm(nrow(dn), 0, 0.3)
  fnl <- suppressWarnings(suppressMessages(
    frm(bf(yn ~ log(c1) * x + exp(a)^k, c1 ~ 1, a ~ 1 + (1 | g), k ~ 1,
           nl = TRUE), family = gaussian(), data = dn,
        prior = set_prior("normal(1.6, 1)", nlpar = "c1"))))
  p <- lap_pair(fnl)
  for (k in c("lap", "full")) {
    p[[k]]$draws[1L, "b_c1_Intercept"] <- -1
    p[[k]]$draws[1L, "b_k_Intercept"] <- 0
  }
  nan <- "NaNs produced"
  full <- allow_warnings(posterior_epred(p$full), nan)
  expect_true(all(is.nan(full[1L, ])))
  expect_true(all(is.finite(full[-1L, ])))
  expect_error(allow_warnings(posterior_epred(p$lap), nan),
               "posterior_epred() needs", fixed = TRUE)
  expect_error(allow_warnings(posterior_linpred(p$lap), nan),
               "posterior_linpred() needs", fixed = TRUE)
  expect_error(allow_warnings(posterior_predict(p$lap), nan),
               "posterior_predict() needs", fixed = TRUE)
  # the guard's condition absent: the same NaN first draw on a quantity
  # the integrated values do not reach is returned, NaN and all
  expect_identical(allow_warnings(posterior_epred(p$lap, dpar = "sigma"),
                                  nan),
                   allow_warnings(posterior_epred(p$full, dpar = "sigma"),
                                  nan))
})

test_that("a group-level slope at a covariate of 0 reads no group value", {
  # Backlog 1530: (0 + x | g) at x = 0 contributes b * 0, which is NA * 0
  # = NA at the probe's NA fill and 0 at every finite one, so the call
  # was refused though it reads no group-level value
  dd <- lap_data()
  fit <- frm(bf(y ~ x + (0 + x | g)), family = gaussian(), data = dd)
  p <- lap_pair(fit)
  nd0 <- data.frame(x = c(0, 0), g = factor(c(1, 2), levels = 1:6))
  for (f in list(posterior_epred, posterior_linpred)) {
    expect_identical(f(p$lap, newdata = nd0), f(p$full, newdata = nd0))
  }
  set.seed(3)
  a <- posterior_predict(p$lap, newdata = nd0)
  set.seed(3)
  expect_identical(a, posterior_predict(p$full, newdata = nd0))
  # a row at x = 0 beside a row that reads the slope is still refused
  nd1 <- data.frame(x = c(0, 1), g = factor(c(1, 2), levels = 1:6))
  expect_error(posterior_epred(p$lap, newdata = nd1),
               "posterior_epred() needs", fixed = TRUE)
  expect_error(posterior_epred(p$lap), "posterior_epred() needs",
               fixed = TRUE)
  # a missing covariate stays NA at every fill, and is not a read
  ndm <- data.frame(x = c(0, NA), g = factor(c(1, 2), levels = 1:6))
  expect_identical(posterior_epred(p$lap, newdata = ndm, re_formula = NA),
                   posterior_epred(p$full, newdata = ndm, re_formula = NA))
})

test_that("a group-level design that cancels or a step() still refuses", {
  # Review 2026-11-11 of lane sampfront, finding 1: a second scalar fill
  # let these through with every effect at 0 (max |laplace - full| 1.01,
  # 1.74, 1.01, 0.30 and 0.93 in dev/review-sampfront-probe.R). Only a
  # design entry that is exactly 0 lets a cell skip the effects.
  set.seed(77)
  G <- 8
  n <- 12
  dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
  dd$x <- stats::rnorm(nrow(dd))
  dd$f <- factor(rep(c("a", "b"), length.out = nrow(dd)))
  u0 <- stats::rnorm(G, 0, 0.8)
  u1 <- stats::rnorm(G, 0, 0.8)
  dd$y <- 0.5 + 0.4 * dd$x + u0[dd$g] + u1[dd$g] * dd$x +
    stats::rnorm(nrow(dd), 0, 0.5)
  q <- function(...) suppressWarnings(suppressMessages(frm(...)))
  msg <- "frm_sample(laplace = TRUE), which integrates them out"
  g1 <- data.frame(g = factor(1, levels = 1:G))
  # A: (1 + x | g) at x = -1, where the two entries are 1 and -1
  pa <- lap_pair(q(bf(y ~ x + (1 + x | g)), family = gaussian(), data = dd))
  nda <- data.frame(x = -1, g = factor(c(1, 2), levels = 1:G))
  expect_error(posterior_epred(pa$lap, newdata = nda), msg, fixed = TRUE)
  expect_error(posterior_linpred(pa$lap, newdata = nda), msg, fixed = TRUE)
  expect_error(posterior_predict(pa$lap, newdata = nda), msg, fixed = TRUE)
  expect_error(conditional_effects(pa$lap, effects = "x", re_formula = NULL,
                                   int_conditions = list(x = -1),
                                   conditions = g1),
               "conditional_effects() needs", fixed = TRUE)
  expect_error(conditional_effects(pa$lap, effects = "x", re_formula = NULL,
                                   method = "posterior_predict",
                                   int_conditions = list(x = -1),
                                   conditions = g1),
               "conditional_effects() needs", fixed = TRUE)
  # B: a sum-coded factor's second level, design row (1, -1)
  # (the design of the group-level term is rebuilt under the session's
  # contrasts, as at fit time, so they stay set for the prediction)
  withr::with_options(list(contrasts = c("contr.sum", "contr.poly")), {
    pb <- lap_pair(q(bf(y ~ f + (1 + f | g)), family = gaussian(),
                     data = dd))
    ndb <- data.frame(f = factor("b", levels = c("a", "b")),
                      g = factor(1:3, levels = 1:G))
    expect_error(posterior_epred(pb$lap, newdata = ndb), msg, fixed = TRUE)
    # the guard's condition absent: the full draws answer this row
    expect_true(all(is.finite(posterior_epred(pb$full, newdata = ndb))))
  })
  # C: (1 + x || g) at x = -1
  pc <- lap_pair(q(bf(y ~ x + (1 + x || g)), family = gaussian(), data = dd))
  expect_error(posterior_epred(pc$lap, newdata = nda), msg, fixed = TRUE)
  # D2: below the population's breakpoint and above the group's, where
  # step() is 0 at fills 0 and 0.5772 alike; a nonlinear body reads its
  # parameters' design only through the body, so it keeps the refusal
  dd$yb <- 1 + 0.8 * pmax(dd$x - (0.2 + u0[dd$g] * 0.5), 0) +
    stats::rnorm(nrow(dd), 0, 0.2)
  fd <- q(bf(yb ~ a + bb * (x - om) * step(x - om), a ~ 1, bb ~ 1,
             om ~ 1 + (1 | g), nl = TRUE), family = gaussian(), data = dd,
          prior = c(set_prior("normal(0, 2)", nlpar = "a"),
                    set_prior("normal(0, 2)", nlpar = "bb"),
                    set_prior("normal(0, 2)", nlpar = "om")))
  pd <- lap_pair(fd)
  beta <- unlist(fd$estimates$beta)
  om0 <- beta[grep("om", names(beta))][1L]
  gmin <- which.min(fd$estimates$b)
  ndd <- data.frame(x = om0 - 0.25, g = factor(gmin, levels = 1:G))
  expect_error(posterior_epred(pd$lap, newdata = ndd), msg, fixed = TRUE)
  # the guard's condition absent: where every entry is 0 the same calls
  # answer, exactly as the full draws do
  pz <- lap_pair(q(bf(y ~ x + (0 + x | g)), family = gaussian(), data = dd))
  nd0 <- data.frame(x = 0, g = factor(c(1, 2), levels = 1:G))
  expect_identical(posterior_epred(pz$lap, newdata = nd0),
                   posterior_epred(pz$full, newdata = nd0))
  ce <- function(d) {
    set.seed(4)
    conditional_effects(d, effects = "x", re_formula = NULL,
                        method = "posterior_predict",
                        int_conditions = list(x = 0), conditions = g1)
  }
  expect_identical(ce(pz$lap)[[1L]]$estimate__, ce(pz$full)[[1L]]$estimate__)
})

## ---- frm_sample() itself --------------------------------------------

test_that("laplace = TRUE with nothing to integrate samples the model", {
  skip_on_cran()
  skip_sampler()
  dd <- lap_data()
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  ess <- c("Effective Samples Size", "R-hat", "Rhat")
  # it used to die in tmbstan on R's "invalid argument to unary operator"
  expect_message(
    a <- allow_warnings(frm_sample(fit, chains = 1, iter = 300,
                                   refresh = 0, seed = 3,
                                   laplace = TRUE), ess),
    "no random effects or mi() values to integrate out", fixed = TRUE)
  b <- allow_warnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3)), ess)
  expect_identical(a$draws, b$draws)
  expect_false(frmtmb.sample:::draws_is_laplace(a))
})

test_that("real laplace draws refuse what needs b and compute the rest", {
  skip_on_cran()
  skip_sampler()
  dd <- lap_data()
  ds <- allow_warnings(suppressMessages(
    frm_sample(bf(y ~ x + ar(t, g) + (1 | g)), family = gaussian(),
               data = dd, chains = 1, iter = 300, refresh = 0, seed = 3,
               laplace = TRUE)),
    c("Effective Samples Size", "R-hat", "Rhat", "lp__"))
  expect_true(frmtmb.sample:::draws_is_laplace(ds))
  expect_error(posterior_predict(ds), "integrates them out", fixed = TRUE)
  ep <- posterior_epred(ds, re_formula = NA)
  expect_true(all(is.finite(ep)))
})

test_that("laplace = TRUE refuses where a REML objective holds beta too", {
  skip_on_cran()
  skip_sampler()
  dd <- lap_data()
  fr <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
            REML = TRUE)
  # sampled as it stands, as check_laplace() samples it, the REML
  # objective integrates the coefficients, and the draws used to carry
  # b_Intercept over theta's values
  expect_error(
    allow_warnings(suppressMessages(
      frm_sample(fr, chains = 1, iter = 200, refresh = 0, seed = 3,
                 laplace = TRUE, .diagnostic = TRUE)), "flat"),
    "would integrate the population-level coefficients", fixed = TRUE)
  # the default priors rebuild the objective with beta sampled, and say
  # that REML's integral over beta is not in it (backlog 1533)
  said <- function(f) {
    ms <- character()
    ds <- allow_warnings(withCallingHandlers(
      frm_sample(f, chains = 1, iter = 200, refresh = 0, seed = 3,
                 laplace = TRUE),
      message = function(m) {
        ms <<- c(ms, conditionMessage(m))
        invokeRestart("muffleMessage")
      }), c("Effective Samples Size", "R-hat", "Rhat", "lp__"))
    list(ds = ds, reml = any(grepl("on a REML model samples the ",
                                   ms, fixed = TRUE)))
  }
  r <- said(fr)
  expect_true(all(c("b_Intercept", "b_x") %in% colnames(r$ds$draws)))
  expect_true(r$reml)
  # the absent case: a maximum-likelihood fit is not told
  expect_false(said(frm(bf(y ~ x + (1 | g)), family = gaussian(),
                        data = dd))$reml)
})

test_that("a one-parameter model samples from its mode", {
  skip_on_cran()
  skip_sampler()
  set.seed(1212)
  dd <- data.frame(x = stats::rnorm(30))
  dd$y <- stats::rpois(30, exp(0.4 * dd$x))
  fx <- frm(bf(y ~ 0 + x), family = poisson(), data = dd)
  expect_length(fx$obj$par, 1L)
  # rstan read a length-one init as a scalar and every chain died with
  # "no more scalars to read"
  ds <- allow_warnings(suppressMessages(
    frm_sample(fx, chains = 2, iter = 300, refresh = 0, seed = 3)),
    c("Effective Samples Size", "R-hat", "Rhat"))
  expect_identical(colnames(ds$draws), c("b_x", "lp__"))
  expect_identical(nrow(ds$draws), 300L)
})
