# re_formula = NA on a model with smooths: what a population-level
# prediction removes and what it keeps.
#
# The rule is brms's: `re_formula` governs the `(... | g)` group-level
# terms and nothing else, so EVERY smooth stays in, a smooth indexed by
# a grouping factor included (`s(g, bs = "re")`, `s(x, g, bs = "fs")`, a
# `t2()` with an `re` margin). Through 0.64.0 those three were dropped
# instead, which answered a different model and made `simulate(NA)` and
# the default `pp_check()` redraw them (dev/resmooth-findings.md).
#
# The reference is mgcv, whose `predict(exclude = )` names one term at a
# time, so the prediction each rule describes can be written down: the
# population one excludes the group-level intercept alone, and the old
# one excluded the factor smooth with it.
#
# Gaps are reported RELATIVE to the size of the values being compared.
# 1e-5 is the level two INDEPENDENT fits agree at, frmtmb's and mgcv's
# (measured: 7.9e-7 on the first block below); 1e-12 is where one fit's
# own two summation orders sit (4.4e-16 absolute on the probe design,
# dev/resmooth-after.txt).

rel_gap <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  max(abs(a - b)) / max(abs(c(a, b)))
}

fosr_data <- function(N = 20L, nt = 15L, seed = 101) {
  set.seed(seed)
  tt <- seq(0, 1, length.out = nt)
  x <- rbinom(N, 1, 0.5)
  b0 <- function(t) 1 + 2 * sin(2 * pi * t)
  bi <- matrix(rnorm(N * 2, 0, 0.6), N, 2)
  Y <- outer(rep(1, N), b0(tt)) + bi[, 1] +
    outer(bi[, 2], tt - 0.5) * 2 +
    matrix(rnorm(N * nt, 0, 0.35), N, nt)
  data.frame(subject = factor(rep(seq_len(N), each = nt)),
             t = rep(tt, N), x = rep(x, each = nt),
             y = as.vector(t(Y)))
}

test_that("re_formula = NA keeps the fs smooth and drops (1 | g)", {
  d <- fosr_data()
  fit <- frm(bf(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5) +
                  (1 | subject)),
             family = gaussian(), data = d)
  gm <- suppressWarnings(
    mgcv::gam(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5) +
                s(subject, bs = "re"), data = d, method = "ML"))

  # conditional prediction: the two packages fit the same model
  expect_lt(max(abs(as.numeric(fitted(fit)[, "Estimate"]) -
                      as.numeric(fitted(gm)))), 1e-4)
  # how far apart two optimizers leave the SAME model on this platform,
  # which is the yardstick for the population comparison below: a fixed
  # 1e-5 passed on Windows and failed at 1.3e-5 on the macOS and Ubuntu
  # CI runners of 2026-09-29, where both fits were right
  cond_gap <- rel_gap(fitted(fit)[, "Estimate"], fitted(gm))

  # population prediction: the group-level intercept goes and both
  # smooths stay, which is mgcv excluding s(subject) alone
  pop <- as.numeric(predict(gm, exclude = "s(subject)"))
  expect_lt(rel_gap(frm_linpred(fit, re_formula = NA), pop),
            10 * max(cond_gap, sqrt(.Machine$double.eps)))

  # SEEN TO FAIL through 0.64.0, which returned this instead: the fs
  # term excluded as well. The two predictions are far apart, so the
  # assertion above cannot pass under the old rule.
  dropped_fs <- as.numeric(predict(gm,
                                   exclude = c("s(t,subject)", "s(subject)")))
  expect_gt(max(abs(pop - dropped_fs)) / stats::sd(pop), 0.1)

  nd <- data.frame(t = seq(0, 1, length.out = 11), x = 0,
                   subject = factor(levels(d$subject)[1],
                                    levels = levels(d$subject)))
  expect_lt(rel_gap(frm_linpred(fit, newdata = nd, re_formula = NA),
                    predict(gm, newdata = nd, exclude = "s(subject)")),
            1e-5)

  # the population prediction still carries a standard error
  se <- frm_linpred(fit, newdata = nd, re_formula = NA, se.fit = TRUE)
  expect_true(all(is.finite(se$se.fit)))
  expect_true(all(se$se.fit > 0))
})

test_that("an fs term needs its grouping column at every re_formula", {
  d <- fosr_data()
  fit <- frm(bf(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5)),
             family = gaussian(), data = d)
  nd <- data.frame(t = seq(0, 1, length.out = 7))   # no `subject` column
  # SEEN TO FAIL through 0.64.0, where re_formula = NA dropped the term
  # and so needed no level. brms refuses the same newdata, and not
  # through its group-level machinery: the grouping factor of a smooth is
  # an ordinary predictor there (dev/resmooth-brms.txt).
  expect_error(frm_linpred(fit, newdata = nd, re_formula = NA),
               "needs the grouping column `subject`")
  expect_error(frm_linpred(fit, newdata = nd),
               "needs the grouping column `subject`")
  # and the refusal no longer offers re_formula = NA as the way out
  expect_error(frm_linpred(fit, newdata = nd),
               "re_formula = NA does not remove the need")

  # a population smooth missing its own covariate is a different fault
  expect_error(frm_linpred(fit, newdata = data.frame(subject = d$subject[1])),
               "needs the column\\(s\\) `t`")

  # with the column present, NA and NULL agree to round-off: the model
  # has no group-level term for NA to drop
  nd2 <- transform(nd, subject = d$subject[1])
  expect_lt(rel_gap(frm_linpred(fit, newdata = nd2, re_formula = NA),
                    frm_linpred(fit, newdata = nd2)), 1e-12)
})

test_that("an unseen fs level errors at every re_formula", {
  d <- fosr_data()
  fit <- frm(bf(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5)),
             family = gaussian(), data = d)
  gm <- suppressWarnings(
    mgcv::gam(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5),
              data = d, method = "ML"))
  nd <- data.frame(t = c(0, 0.5, 1), subject = factor("brand_new"))
  expect_error(frm_linpred(fit, newdata = nd),
               "New levels in the factor-smooth term")
  # SEEN TO FAIL through 0.64.0: re_formula = NA answered here
  expect_error(frm_linpred(fit, newdata = nd, re_formula = NA),
               "New levels in the factor-smooth term")
  # allow_new_levels is the only way out, and it gives the population
  # curve because mgcv's fs basis zero-rows a level it does not know
  expect_lt(rel_gap(frm_linpred(fit, newdata = nd, allow_new_levels = TRUE),
                    predict(gm,
                            newdata = transform(nd,
                                                subject = d$subject[1]),
                            exclude = "s(t,subject)")),
            1e-5)
})

test_that("s(g, bs = 're') is a smooth, so re_formula = NA keeps it", {
  set.seed(7)
  n <- 240L
  d <- data.frame(t = runif(n), g = factor(sample(letters[1:8], n, TRUE)))
  d$y <- sin(3 * d$t) + rnorm(8, 0, 0.7)[d$g] + rnorm(n, 0, 0.3)

  fre <- frm(bf(y ~ s(t, k = 6) + s(g, bs = "re")),
             family = gaussian(), data = d)
  gre <- mgcv::gam(y ~ s(t, k = 6) + s(g, bs = "re"), data = d,
                   method = "ML")
  # SEEN TO FAIL through 0.64.0, which excluded s(g): the two blocks of
  # the term are its wiggly part and its one unpenalized column, and the
  # old rule kept block 1 and dropped block 2
  expect_lt(rel_gap(frm_linpred(fre, re_formula = NA), predict(gre)), 1e-5)
  expect_lt(rel_gap(frm_linpred(fre, re_formula = NA),
                    frm_linpred(fre)), 1e-12)
  # the whole term, not just its unpenalized part: the old answer
  # differs by more than the group spread
  pop_old <- as.numeric(predict(gre, exclude = "s(g)"))
  expect_gt(max(abs(as.numeric(predict(gre)) - pop_old)) / stats::sd(d$y),
            0.1)
  # the wiggly part of s(t) survived: a flat line would not
  expect_gt(diff(range(frm_linpred(fre, re_formula = NA))), 0.5)

  # an unseen level of a factor bs = "re" smooth: named refusal by
  # default, and a second named refusal under allow_new_levels = TRUE
  # (the design has one column per fitted level, no zero row).
  # re_formula = NA is no longer a way out of either, because it keeps
  # the term: SEEN TO FAIL through 0.64.0, where it answered with 2
  # values.
  nd_new <- data.frame(t = c(0.2, 0.6), g = factor(c("zz", "a")))
  expect_error(frm_linpred(fre, newdata = nd_new), "New levels")
  expect_error(frm_linpred(fre, newdata = nd_new, allow_new_levels = TRUE),
               "no zero row")
  expect_error(frm_linpred(fre, newdata = nd_new, re_formula = NA),
               "New levels")
})

test_that("a t2() with an re margin is kept by re_formula = NA", {
  set.seed(11)
  n <- 240L
  d <- data.frame(t = runif(n), g = factor(sample(letters[1:8], n, TRUE)))
  d$y <- sin(3 * d$t) + rnorm(8, 0, 0.7)[d$g] * d$t + rnorm(n, 0, 0.3)
  fit <- suppressWarnings(frm(bf(y ~ t2(t, g, bs = c("cr", "re"), k = 5)),
                              family = gaussian(), data = d))
  gm <- suppressWarnings(mgcv::gam(y ~ t2(t, g, bs = c("cr", "re"), k = 5),
                                   data = d, method = "ML"))
  # every block of the term is kept, so NA is the full prediction.
  # SEEN TO FAIL through 0.64.0, which dropped both blocks and left the
  # intercept alone.
  # 1e-5, the level the file header states for two independent fits;
  # measured 2.277e-08 here (dev/resmooth-t2gap.R)
  expect_lt(rel_gap(frm_linpred(fit, re_formula = NA), predict(gm)), 1e-5)
  expect_lt(rel_gap(frm_linpred(fit, re_formula = NA), frm_linpred(fit)),
            1e-12)
  # the old answer was the intercept at every row
  expect_gt(diff(range(frm_linpred(fit, re_formula = NA))) / stats::sd(d$y),
            0.5)
})

test_that("the group/population split is read off the smooth object", {
  set.seed(9)
  d <- data.frame(t = runif(120), x = rnorm(120),
                  g = factor(sample(letters[1:5], 120, TRUE)))
  cl <- function(spec) {
    mgcv::smoothCon(spec, data = d, absorb.cons = TRUE, modCon = 3)[[1L]]
  }
  gv <- function(spec) frmtmb:::smooth_group_var(cl(spec), d)
  expect_identical(gv(mgcv::s(t, g, bs = "fs", k = 4)), "g")
  expect_identical(gv(mgcv::s(g, bs = "re")), "g")
  expect_identical(gv(mgcv::s(x, g, bs = "re")), "g")
  expect_identical(gv(mgcv::t2(t, g, bs = c("cr", "re"), k = 4)), "g")
  expect_null(gv(mgcv::s(t, k = 5)))
  expect_null(gv(mgcv::s(t, by = x, k = 5)))
  expect_null(gv(mgcv::s(t, by = g, k = 5)))
  # sz writes its level curves as contrasts against a reference level,
  # so it is a fixed effect even though it names a factor the way fs does
  expect_null(gv(mgcv::s(t, g, bs = "sz", k = 4)))
})

test_that("conditional_effects draws the reference level's fs curve", {
  d <- fosr_data()
  fit <- frm(bf(y ~ x + s(t, k = 8) + s(t, subject, bs = "fs", k = 5)),
             family = gaussian(), data = d)
  gm <- suppressWarnings(
    mgcv::gam(y ~ x + s(t, k = 8) + s(t, subject, bs = "fs", k = 5),
              data = d, method = "ML"))
  ce <- conditional_effects(fit, effects = "t", resolution = 9)
  df <- ce[["t"]]
  # the display holds the other predictors at their reference values,
  # `subject` among them, and its curve is now part of the answer:
  # conditional_effects() passes re_formula = NA, which keeps every
  # smooth. SEEN TO FAIL through 0.64.0, where the term was dropped and
  # this fit's only other term was the intercept, so the drawn curve was
  # FLAT (range 0.4901 to 0.4901 on the probe design,
  # dev/resmooth-before.txt).
  nd <- data.frame(t = df$t, x = mean(d$x),
                   subject = factor(levels(d$subject)[1],
                                    levels = levels(d$subject)))
  ref <- as.numeric(predict(gm, newdata = nd))
  expect_lt(rel_gap(df$estimate__, ref), 1e-5)

  # the grouping factor is not offered as an effect to draw
  expect_false("subject" %in% names(conditional_effects(fit)))
  expect_true("t" %in% names(conditional_effects(fit)))
})

test_that("conditional_effects names matrix columns as the reason", {
  set.seed(202)
  n <- 120L
  nS <- 20L
  S <- seq(0, 1, length.out = nS)
  Xf <- t(replicate(n, cumsum(rnorm(nS, 0, 0.4)) + rnorm(nS, 0, 0.2)))
  sof <- data.frame(y = 1 + as.vector((Xf / nS) %*% (2 * sin(2 * pi * S))) +
                      rnorm(n, 0, 0.4))
  sof$Smat <- matrix(S, n, nS, byrow = TRUE)
  sof$LX <- Xf / nS
  fit <- frm(bf(y ~ s(Smat, by = LX, k = 8)), family = gaussian(),
             data = sof)
  expect_error(conditional_effects(fit), "matrix column\\(s\\)")
  expect_error(conditional_effects(fit), "`Smat`")
  expect_error(conditional_effects(fit), "predict\\(newdata = \\)")
  expect_error(conditional_effects(fit), "s\\(Smat\\):LX")

  # the generic refusal is still the one an empty predictor gets
  set.seed(3)
  d0 <- data.frame(y = rnorm(60))
  fit0 <- frm(bf(y ~ 1), family = gaussian(), data = d0)
  expect_error(conditional_effects(fit0), "No plottable predictors found")
})
