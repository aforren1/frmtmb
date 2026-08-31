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
