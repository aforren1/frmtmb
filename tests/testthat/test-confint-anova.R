fit_ca <- local({
  data(sleepstudy, package = "lme4")
  list(
    data = sleepstudy,
    full = frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
                  data = sleepstudy),
    null = frm(bf(Reaction ~ 1 + (Days | Subject)) + gaussian(),
                  data = sleepstudy)
  )
})

test_that("wald confidence intervals are consistent with vcov", {
  fit <- fit_ca$full
  ci <- confint(fit)
  expect_true(all(c("(Intercept)", "Days", "sigma_(Intercept)") %in%
                    rownames(ci)))
  expect_true(all(ci[, "lwr"] < ci[, "est"] & ci[, "est"] < ci[, "upr"]))
  se_days <- sqrt(vcov(fit)["Days", "Days"])
  expect_lt(abs((ci["Days", "upr"] - ci["Days", "lwr"]) / 2 -
                  qnorm(0.975) * se_days), 1e-8)
})

test_that("profile and uniroot intervals agree with each other", {
  fit <- fit_ca$full
  cp <- confint(fit, parm = "Days", method = "profile")
  cu <- confint(fit, parm = "Days", method = "uniroot")
  expect_lt(abs(cp["Days", "lwr"] - cu["Days", "lwr"]), 1e-2)
  expect_lt(abs(cp["Days", "upr"] - cu["Days", "upr"]), 1e-2)
  # profile interval is close to wald for this near-quadratic likelihood
  cw <- confint(fit, parm = "Days")
  expect_lt(abs(cp["Days", "lwr"] - cw["Days", "lwr"]), 0.2)
})

test_that("anova performs likelihood-ratio tests", {
  tab <- anova(fit_ca$null, fit_ca$full)
  expect_s3_class(tab, "anova")
  expect_identical(nrow(tab), 2L)
  expect_gt(tab$Chisq[2], 0)
  expect_identical(tab[["Chi Df"]][2], 1L)
  expect_lt(tab[["Pr(>Chisq)"]][2], 1e-4)

  reml <- update(fit_ca$full, REML = TRUE)
  expect_error(anova(reml, fit_ca$null), "REML")
})

test_that("anova(refit = TRUE) compares REML fits as ML fits", {
  r_full <- update(fit_ca$full, REML = TRUE)
  r_null <- update(fit_ca$null, REML = TRUE)
  # the fixed designs differ, so the two restricted likelihoods are for
  # different error contrasts and the REML path refuses
  expect_error(anova(r_full, r_null), "same column space")
  expect_message(tab <- anova(r_full, r_null, refit = TRUE),
                 "refitting 2 REML models with ML")
  ref <- anova(fit_ca$full, fit_ca$null)
  expect_identical(rownames(tab), rownames(ref))
  expect_identical(tab$Df, ref$Df)
  expect_equal(tab$logLik, ref$logLik, tolerance = 1e-6)
  expect_equal(tab$Chisq, ref$Chisq, tolerance = 1e-6)
  # a REML/ML mix is refused by default and refit under refit = TRUE
  expect_error(anova(r_full, fit_ca$null), "cannot mix REML and ML")
  expect_message(anova(r_full, fit_ca$null, refit = TRUE),
                 "refitting 1 REML model with ML")
})

test_that("update re-fits with modified arguments", {
  reml <- update(fit_ca$full, REML = TRUE)
  expect_true(reml$REML)
  expect_false(isTRUE(all.equal(as.numeric(logLik(reml)),
                                as.numeric(logLik(fit_ca$full)))))
  cl <- update(fit_ca$full, REML = TRUE, evaluate = FALSE)
  expect_true(is.call(cl))
})

test_that("diagnose reports a clean fit as clean", {
  d <- diagnose(fit_ca$full, quiet = TRUE)
  expect_equal(d$convergence, 0)
  expect_true(d$pdHess)
  expect_lt(d$max_grad, 1e-2)
  expect_length(d$bad_se, 0)
  expect_output(diagnose(fit_ca$full), "No convergence problems")
})

# --- anova() and nesting (the habit-model finding) --------------------

test_that("anova() warns when the fixed effects are not nested", {
  set.seed(77)
  dd <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150))
  dd$y <- stats::rnorm(150, 1 + 0.6 * dd$x, 1)
  mx <- frm(bf(y ~ x) + gaussian(), data = dd)
  mz <- frm(bf(y ~ z) + gaussian(), data = dd)
  expect_warning(anova(mx, mz), "are not nested")
  # the table is still produced: the warning is a caveat, not a refusal
  tab <- suppressWarnings(anova(mx, mz))
  expect_s3_class(tab, "anova")
})

test_that("anova() stays silent on nested pairs, including a rebasis", {
  set.seed(78)
  dd <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150),
                   g = factor(rep(1:15, 10)))
  dd$y <- stats::rnorm(150, 1 + 0.6 * dd$x +
                         stats::rnorm(15, 0, 0.5)[dd$g], 1)
  m0 <- frm(bf(y ~ x) + gaussian(), data = dd)
  m1 <- frm(bf(y ~ x + z) + gaussian(), data = dd)
  expect_no_warning(anova(m0, m1))
  # a variance component added: same fixed effects, so nested
  m2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_no_warning(anova(m0, m2))
  # different names for the same span, plus a term: the column-space
  # fallback is what keeps this quiet
  mp <- frm(bf(y ~ poly(x, 2)) + gaussian(), data = dd)
  expect_no_warning(anova(m0, mp))
})

test_that("a dpar predictor counts as a fixed effect for nesting", {
  set.seed(79)
  dd <- data.frame(x = stats::rnorm(200), z = stats::rnorm(200),
                   w = stats::rnorm(200))
  dd$y <- stats::rnorm(200, 1 + 0.5 * dd$x, exp(0.2 + 0.3 * dd$z))
  a <- frm(bf(y ~ x, sigma ~ z) + gaussian(), data = dd)
  b <- frm(bf(y ~ x, sigma ~ w) + gaussian(), data = dd)
  expect_warning(anova(a, b), "are not nested")
  c0 <- frm(bf(y ~ x) + gaussian(), data = dd)
  expect_no_warning(anova(c0, a))
})
