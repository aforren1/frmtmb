# conditional_effects() and emmeans() on a model with an offset().
# Before, the model frame held the column `offset(log(time))` and not
# `time`, so neither grid had `time`: log(time) then found stats::time()
# ("non-numeric argument to mathematical function") and emmeans()'s
# grid stopped with "undefined columns selected" (dev/formrobust-repro2.R
# for the before and after arms). brms holds `time` at its mean in both
# grids (prepare_conditions()). Its emmeans() adds a predictor's offset
# back through emmeans's own `.offset.` column at the grid's value, and
# leaves out an offset inside a nonlinear parameter's formula, which its
# reference grid cannot see (the review's brms runs,
# dev/formrobust-rev-log/brms-emm.txt and brms-emm2.txt, and the block
# at fixed parameters at the end of this file).

og_data <- local({
  set.seed(21)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), time = stats::runif(n, 1, 3),
                  f = gl(2, 40), z = stats::runif(n))
  d$yc <- stats::rpois(n, exp(0.2 + 0.3 * d$x) * d$time)
  d$y <- stats::rnorm(n, 0.5 * d$x + d$z, 1)
  d$ys <- stats::rnorm(n, 0.5 * d$x, exp(0.2 + log(d$time)))
  d
})

og_rel <- function(a, b) max(abs(a - b) / pmax(abs(a), abs(b)))

test_that("conditional_effects() holds an offset's variable at its mean", {
  fit <- frm(bf(yc ~ x + f + offset(log(time))), data = og_data,
             family = poisson())
  b <- fixef(fit)[, "Estimate"]
  ce <- conditional_effects(fit, effects = "x")[[1L]]
  ref <- exp(b[["Intercept"]] + b[["x"]] * ce$x + log(mean(og_data$time)))
  # identity up to rounding in the exp(): judged in ulps of the value
  expect_lt(og_rel(ce$estimate__, ref), 64 * .Machine$double.eps)
  # an offset is not a default display, as brms's get_all_effects()
  expect_identical(names(conditional_effects(fit)), c("x", "f"))
})

test_that("an offset's variable is not a valid effect, as in brms", {
  # brms:::get_all_effects() on this formula is x and f alone, so brms
  # stops on `time` and drops it beside a valid effect. The frame holds
  # `time` for the grid, which let it through before the 0.67.0
  # consolidation (dev/rel067-effoffset.R, dev/rel067-log/effoffset.txt)
  fit <- frm(bf(yc ~ x + f + offset(log(time))), data = og_data,
             family = poisson())
  expect_error(conditional_effects(fit, effects = "time"),
               "All specified effects are invalid for this model")
  ce <- allow_warnings(conditional_effects(fit, effects = c("x", "time")),
                       "Some specified effects are invalid",
                       require = "Some specified effects are invalid")
  expect_identical(names(ce), "x")
  # the guard-absent case: a variable the offset and a term both read
  fb <- frm(bf(yc ~ x + time + offset(log(time))), data = og_data,
            family = poisson())
  expect_identical(names(conditional_effects(fb, effects = "time")),
                   "time")
})

test_that("emmeans() includes the offset at the grid's value, as brms", {
  skip_if_not_installed("emmeans")
  fit <- frm(bf(yc ~ x + f + offset(log(time))), data = og_data,
             family = poisson())
  b <- fixef(fit)[, "Estimate"]
  lin <- b[["Intercept"]] + b[["x"]] * mean(og_data$x) + c(0, b[["f2"]])
  lmt <- log(mean(og_data$time))
  em <- summary(emmeans::emmeans(fit, ~ f))
  expect_lt(og_rel(em$emmean, lin + lmt), 64 * .Machine$double.eps)
  # the guard-absent case: at time = 1 the offset is log(1) = 0
  e1 <- summary(emmeans::emmeans(fit, ~ f, at = list(time = 1)))
  expect_lt(og_rel(e1$emmean, lin), 64 * .Machine$double.eps)
  # epred = TRUE includes the offset once
  ep <- summary(emmeans::emmeans(fit, ~ f, epred = TRUE))
  expect_lt(og_rel(ep$emmean, exp(lin + lmt)), 64 * .Machine$double.eps)
})

test_that("a dpar's offset is included, a nonlinear parameter's is not", {
  skip_if_not_installed("emmeans")
  lmt <- log(mean(og_data$time))
  fs <- frm(bf(ys ~ x, sigma ~ 1 + offset(log(time))), data = og_data)
  es <- summary(emmeans::emmeans(fs, ~ 1, dpar = "sigma"))
  expect_lt(og_rel(es$emmean, fixef(fs)["sigma_Intercept", 1] + lmt),
            64 * .Machine$double.eps)
  # an offset inside `a ~ ...` of a nonlinear body is left out, as brms's
  # grid leaves it out; mean(z) is 0.5 here, so keeping it would show
  fn <- frm(bf(y ~ a + b * x, a ~ 1 + offset(z), b ~ 1, nl = TRUE),
            data = og_data)
  bn <- fixef(fn)[, "Estimate"]
  en <- summary(emmeans::emmeans(fn, ~ 1))
  expect_lt(og_rel(en$emmean, bn[["a_Intercept"]] +
                     bn[["b_Intercept"]] * mean(og_data$x)),
            64 * .Machine$double.eps)
  # the parameter itself, selected with nlpar =, keeps its offset
  ea <- summary(emmeans::emmeans(fn, ~ 1, nlpar = "a"))
  expect_lt(og_rel(ea$emmean, bn[["a_Intercept"]] + mean(og_data$z)),
            64 * .Machine$double.eps)
})

test_that("the grid route and a multivariate fit add each offset once", {
  skip_if_not_installed("emmeans")
  lmt <- log(mean(og_data$time))
  # re_formula = NULL takes the grid route, whose prediction carries the
  # offset itself; emmeans's .offset. must not add it a second time
  d <- og_data
  d$g <- gl(8, 10)
  fr <- frm(bf(yc ~ x + (1 | g) + offset(log(time))), data = d,
            family = poisson())
  eg <- summary(emmeans::emmeans(fr, ~ g, re_formula = NULL))
  nd <- data.frame(x = mean(d$x), time = mean(d$time), g = levels(d$g))
  expect_lt(og_rel(eg$emmean, as.vector(frm_linpred(fr, newdata = nd,
                                                    type = "link"))),
            64 * .Machine$double.eps)
  # both responses of a multivariate fit on one grid: each its own
  # offset, where emmeans's single .offset. column would give both
  fm <- frm(bf(y ~ x + offset(z)) + bf(ys ~ x), data = og_data)
  bm <- fixef(fm)[, "Estimate"]
  em <- summary(emmeans::emmeans(fm, ~ rep.meas))
  mx <- mean(og_data$x)
  expect_lt(og_rel(em$emmean,
                   c(bm[["y_Intercept"]] + bm[["y_x"]] * mx +
                       mean(og_data$z),
                     bm[["ys_Intercept"]] + bm[["ys_x"]] * mx)),
            64 * .Machine$double.eps)
})

test_that("an offset in other positions reaches both grids", {
  skip_if_not_installed("emmeans")
  fits <- list(
    frm(bf(yc ~ x + f + offset(time)), data = og_data, family = poisson()),
    frm(bf(yc ~ x + f + offset(log(time) + 0.1)), data = og_data,
        family = poisson()),
    frm(bf(y ~ x + f, sigma ~ x + offset(log(time))), data = og_data),
    frm(bf(y ~ a + b * x, a ~ f + offset(z), b ~ 1, nl = TRUE),
        data = og_data)
  )
  for (fit in fits) {
    ce <- conditional_effects(fit, effects = "x")[[1L]]
    expect_true(all(is.finite(ce$estimate__)))
    em <- summary(emmeans::emmeans(fit, ~ f))
    expect_true(all(is.finite(em$emmean)))
  }
  # a written constant in the offset moves the intercept by as much, and
  # the emmean, which includes the offset, does not move: the same model
  f0 <- frm(bf(yc ~ x + f + offset(log(time))), data = og_data,
            family = poisson())
  s0 <- summary(emmeans::emmeans(f0, ~ f))
  e1 <- summary(emmeans::emmeans(fits[[2L]], ~ f))$emmean
  expect_lt(max(abs(e1 - s0$emmean) / s0$SE), 1e-4)
  # a factor read inside an offset is not a factor of the terms, and
  # model.frame() warned that it "is not a factor"
  d <- og_data
  d$ef <- factor(ifelse(d$time > 2, "hi", "lo"))
  ff <- frm(bf(yc ~ x + offset(log(as.numeric(ef)))), data = d,
            family = poisson())
  expect_no_warning(conditional_effects(ff, effects = "x"))
  expect_no_warning(emmeans::emmeans(ff, ~ 1, epred = TRUE))
})

test_that("the model frame holds an offset's variable", {
  fit <- frm(bf(yc ~ x + offset(log(time))), data = og_data,
             family = poisson())
  expect_true(all(c("offset(log(time))", "time") %in%
                    names(model.frame(fit))))
  # the guard-absent case: an offset whose variable is a constant of the
  # formula environment adds no column, and the fit is the offset's
  k <- 0.5
  f1 <- frm(bf(yc ~ x + offset(log(time) * k)), data = og_data,
            family = poisson())
  expect_false("k" %in% names(model.frame(f1)))
})

test_that("emmeans() equals brms's emmeans() at fixed parameters", {
  skip_unless_brms_fit()
  skip_if_not_installed("emmeans")
  # brms's draws are frmtmb's estimates (Fixed_param), so brms's
  # emmean, a median over identical draws, is the same number
  cmp <- function(bform, fit, family, specs, ...) {
    bb <- brms_fixed_fit(bform, family, og_data, fit, ndraws = 4)
    a <- summary(emmeans::emmeans(bb, specs, ...))$emmean
    b <- summary(emmeans::emmeans(fit, specs, ...))$emmean
    expect_exact_num(a, b, label = deparse1(bform$formula))
  }
  fp <- frm(bf(yc ~ x + f + offset(log(time))), data = og_data,
            family = poisson())
  cmp(brms::bf(yc ~ x + f + offset(log(time))), fp, poisson(), ~ f)
  cmp(brms::bf(yc ~ x + f + offset(log(time))), fp, poisson(), ~ f,
      at = list(time = 1))
  cmp(brms::bf(yc ~ x + f + offset(log(time))), fp, poisson(), ~ f,
      epred = TRUE)
  fs <- frm(bf(ys ~ x, sigma ~ 1 + offset(log(time))), data = og_data)
  cmp(brms::bf(ys ~ x, sigma ~ 1 + offset(log(time))), fs, gaussian(),
      ~ 1, dpar = "sigma")
})
