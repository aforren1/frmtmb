# Two defects of item 2.6d/2.6f found at the 0.61.0 consolidation, both
# on fits whose fixed effects are integrated out (REML = TRUE, or
# control(profile = TRUE)): there the fixed effects are not outer
# parameters, and code that read the OUTER covariance by position, or
# drew over the outer vector alone, went wrong. The ML paths are
# unchanged, bitwise.

reml_fixture <- function() {
  set.seed(1)
  d <- data.frame(x = stats::rnorm(40), g = factor(rep(1:8, 5)))
  d$y <- 1 + d$x + stats::rnorm(8)[d$g] + stats::rnorm(40)
  d
}

test_that("summary() answers on REML and profile fits", {
  d <- reml_fixture()
  for (f in list(frm(y ~ x + (1 | g), data = d, REML = TRUE),
                 frm(y ~ x + (1 | g), data = d,
                     control = frmtmb_control(profile = TRUE)))) {
    # it stopped with "row names contain missing values": sigma's
    # standard error was read from a covariance with no fixed-effect rows
    s <- summary(f)
    expect_true(is.finite(s$spec_pars["sigma", "Est.Error"]))
    # and it is the delta-method standard error of exp(log sigma), read
    # from the covariance over the estimated coefficients
    cf <- frmtmb:::fixef_estimated(f)
    V <- vcov_estimated(f)
    k <- length(cf)
    expect_equal(s$spec_pars["sigma", "Est.Error"],
                 exp(cf[[k]]) * sqrt(V[k, k]), tolerance = 1e-6)
  }
})

test_that("predict() draws the fixed effects under REML", {
  d <- reml_fixture()
  f <- frm(y ~ x + (1 | g), data = d, REML = TRUE)
  # the draw space carries beta, with the covariance vcov_estimated()
  # reports for it
  ds <- frmtmb:::fit_draw_space(f)
  nb <- sum(ds$map$comp == "beta")
  expect_equal(nb, 2L)
  expect_equal(unname(ds$V[seq_len(nb), seq_len(nb)]),
               unname(vcov_estimated(f)[1:2, 1:2]), tolerance = 1e-10)
  # at an extrapolated point the fixed-effect uncertainty is most of the
  # interval: the draws left it out and gave 1.0147 against the analytic
  # sqrt(se^2 + sigma^2) = 1.1611, 13% narrow
  nd <- data.frame(x = 3)
  se_lin <- frm_linpred(f, newdata = nd, re_formula = NA,
                        se.fit = TRUE)$se.fit
  sig <- exp(frmtmb:::fixef_estimated(f)[[3]])
  set.seed(9)
  p <- predict(f, newdata = nd, re_formula = NA, ndraws = 4000)
  ratio <- p[1, "Est.Error"] / sqrt(se_lin^2 + sig^2)
  expect_gt(ratio, 0.97)
  expect_lt(ratio, 1.06)
})
