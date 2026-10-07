# A prediction along a lost direction on a fit WITH random effects.
#
# Lane nanse made a prediction that moves along a direction the fit
# does not determine report NaN and one warning, but only on fits
# without random effects: with them, predictions read the joint
# precision, which carried no lost directions. On this fixture (the
# OLRE fit of test-se-check.R: one observation per id, so the id sd and
# sigma are confounded and both lose their standard errors) 0.68.0
# reported sigma's standard error as 0 0 0 with no warning, under ML and
# under REML. Since 0.68.1 get_joint_cov() carries the lost directions
# with random effects too (joint_cov_repair(); dev/cifix-findings.md).

olre_fit <- function(reml) {
  set.seed(5)
  d <- data.frame(id = factor(1:80), x = rnorm(80))
  d$y <- 1 + 0.5 * d$x + rnorm(80)
  fit <- allow_warnings(
    frm(bf(y ~ x + (1 | id)), family = gaussian(), data = d, REML = reml,
        control = frmtmb_control(check_olre = "ignore")),
    "Standard errors are not available")
  list(fit = fit, nd = d[1:3, ])
}

lost_warnings <- function(expr) {
  n <- 0L
  val <- withCallingHandlers(expr, frmtmb_se_lost_prediction = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  list(value = val, n = n)
}

test_that("sigma along a lost direction is NaN with one warning, ML and REML", {
  for (reml in c(FALSE, TRUE)) {
    o <- olre_fit(reml)
    lab <- paste("REML =", reml)
    lost <- names(sdr_of(o$fit)$se_lost)
    # the premise: sigma's coefficient is among the lost parameters
    expect_true("sigma_(Intercept)" %in% lost, info = lab)
    s <- lost_warnings(frm_linpred(o$fit, newdata = o$nd, dpar = "sigma",
                                   se.fit = TRUE, re_formula = NA))
    expect_true(all(is.nan(s$value$se.fit)), info = lab)
    expect_identical(s$n, 1L, info = lab)
    f <- lost_warnings(fitted(o$fit, newdata = o$nd, dpar = "sigma",
                              re_formula = NA))
    expect_true(all(is.nan(f$value[, "Est.Error"])), info = lab)
    expect_identical(f$n, 1L, info = lab)
    # mu does not load sigma's coefficient and keeps its standard error,
    # without a warning
    m <- lost_warnings(frm_linpred(o$fit, newdata = o$nd, dpar = "mu",
                                   se.fit = TRUE, re_formula = NA))
    expect_true(all(is.finite(m$value$se.fit) & m$value$se.fit > 0),
                info = lab)
    expect_identical(m$n, 0L, info = lab)
  }
})

test_that("frm_joint_cov() shows the lost parameters as NaN", {
  # Moved here from test-se-check.R, whose smooth fit loses its sds
  # with some BLAS builds and not others (dev/ciharden-findings.md);
  # this fit loses two by construction, on every platform
  o <- olre_fit(FALSE)
  lost <- sdr_of(o$fit)$se_lost
  expect_length(lost, 2L)
  jc <- frm_joint_cov(o$fit)
  bad <- jc$lost_pos
  expect_length(bad, length(lost))
  # one a fixed coefficient (sigma's) and one an outer parameter (the
  # id sd), so the positions are read from the joint layout, not
  # matched by name
  expect_setequal(jc$names[bad], c("betad", "theta"))
  expect_true(all(is.nan(jc$V[bad, ])))
  expect_true(all(is.nan(jc$V[, bad])))
  expect_true(all(is.finite(jc$V[-bad, -bad])))
})
