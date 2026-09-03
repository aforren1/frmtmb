# re.form = NA on a model with smooths: what is a population effect and
# what is a group deviation.
#
# The reference is mgcv, which spells the same distinction as
# predict(gam, exclude = ). A smooth's wiggly part is a random-effect
# block in the fitted objective, so the two families of block have to be
# told apart by what they MEAN, not by where their coefficients live.

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

test_that("re.form = NA keeps the population smooth and drops the fs term", {
  d <- fosr_data()
  fit <- frm(bf(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5) +
                  (1 | subject)),
             family = gaussian(), data = d)
  gm <- suppressWarnings(
    mgcv::gam(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5) +
                s(subject, bs = "re"), data = d, method = "ML"))

  # conditional prediction: the two packages fit the same model
  expect_lt(max(abs(as.numeric(fitted(fit)) -
                      as.numeric(fitted(gm)))), 1e-4)

  # population prediction: mgcv's own exclusion of the two group terms
  pop_gam <- as.numeric(predict(gm,
                                exclude = c("s(t,subject)", "s(subject)")))
  expect_lt(max(abs(as.numeric(predict(fit, re.form = NA)) - pop_gam)), 1e-6)

  # and it is NOT merely the fs term left in (the pre-fix behavior)
  kept_fs <- as.numeric(predict(gm, exclude = "s(subject)"))
  expect_gt(max(abs(pop_gam - kept_fs)), 0.1)

  nd <- data.frame(t = seq(0, 1, length.out = 11), x = 0,
                   subject = factor(levels(d$subject)[1],
                                    levels = levels(d$subject)))
  expect_lt(max(abs(
    as.numeric(predict(fit, newdata = nd, re.form = NA)) -
      as.numeric(predict(gm, newdata = nd,
                         exclude = c("s(t,subject)", "s(subject)"))))), 1e-6)

  # the population prediction still carries a standard error
  se <- predict(fit, newdata = nd, re.form = NA, se.fit = TRUE)
  expect_true(all(is.finite(se$se.fit)))
  expect_true(all(se$se.fit > 0))
})

test_that("re.form = NA needs no grouping column for a dropped fs term", {
  d <- fosr_data()
  fit <- frm(bf(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5)),
             family = gaussian(), data = d)
  gm <- suppressWarnings(
    mgcv::gam(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5),
              data = d, method = "ML"))
  nd <- data.frame(t = seq(0, 1, length.out = 7))   # no `subject` column
  expect_lt(max(abs(
    as.numeric(predict(fit, newdata = nd, re.form = NA)) -
      as.numeric(predict(gm,
                         newdata = transform(nd,
                                             subject = d$subject[1]),
                         exclude = "s(t,subject)")))), 1e-6)

  # the conditional prediction does need it, and says so by name
  expect_error(predict(fit, newdata = nd),
               "needs the grouping column `subject`")
  expect_error(predict(fit, newdata = nd), "re\\.form = NA")

  # a population smooth missing its own covariate is a different fault
  expect_error(predict(fit, newdata = data.frame(subject = d$subject[1])),
               "needs the column\\(s\\) `t`")
})

test_that("an unseen fs level errors, and is allowed at the population level", {
  d <- fosr_data()
  fit <- frm(bf(y ~ s(t, k = 8) + s(t, subject, bs = "fs", k = 5)),
             family = gaussian(), data = d)
  nd <- data.frame(t = c(0, 0.5, 1), subject = factor("brand_new"))
  expect_error(predict(fit, newdata = nd),
               "New levels in the factor-smooth term")
  # allowed: the term contributes nothing, which is the population curve
  expect_equal(as.numeric(predict(fit, newdata = nd,
                                  allow_new_levels = TRUE)),
               as.numeric(predict(fit, newdata = nd, re.form = NA)),
               tolerance = 1e-10)
})

test_that("s(g, bs = 're') is a group-level smooth", {
  set.seed(7)
  n <- 240L
  d <- data.frame(t = runif(n), g = factor(sample(letters[1:8], n, TRUE)))
  d$y <- sin(3 * d$t) + rnorm(8, 0, 0.7)[d$g] + rnorm(n, 0, 0.3)

  fre <- frm(bf(y ~ s(t, k = 6) + s(g, bs = "re")),
             family = gaussian(), data = d)
  gre <- mgcv::gam(y ~ s(t, k = 6) + s(g, bs = "re"), data = d,
                   method = "ML")
  expect_lt(max(abs(as.numeric(predict(fre, re.form = NA)) -
                      as.numeric(predict(gre, exclude = "s(g)")))), 1e-5)
  # the wiggly part of s(t) survived: a flat line would not
  expect_gt(diff(range(predict(fre, re.form = NA))), 0.5)
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

test_that("conditional_effects draws the population smooth on an fs fit", {
  d <- fosr_data()
  fit <- frm(bf(y ~ x + s(t, k = 8) + s(t, subject, bs = "fs", k = 5)),
             family = gaussian(), data = d)
  gm <- suppressWarnings(
    mgcv::gam(y ~ x + s(t, k = 8) + s(t, subject, bs = "fs", k = 5),
              data = d, method = "ML"))
  ce <- conditional_effects(fit, effects = "t", resolution = 9)
  df <- ce[["t"]]
  # the display holds the other predictors at their reference values
  nd <- data.frame(t = df$t, x = mean(d$x),
                   subject = factor(levels(d$subject)[1],
                                    levels = levels(d$subject)))
  ref <- as.numeric(predict(gm, newdata = nd, exclude = "s(t,subject)"))
  expect_lt(max(abs(df$estimate__ - ref)), 1e-5)

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
