# A written theta formula's reference component, as brms 2.23.0 sets it
# (dev/correct-log/brms-theta.txt, script dev/correct-brms-probe.R):
# brms predicts the mixing weights only when all but one theta has a
# formula, and the unwritten one is the reference, held at 0:
#   K = 2, theta2 ~ x:        vector[N] theta1 = rep_vector(0.0, N)
#   K = 3, theta2, theta3 ~ x: theta1 = rep_vector(0.0, N)
#   K = 3, theta1 ~ x:        "Can only predict all but one mixing
#                             proportion."
# Through 0.61.0 frmtmb held the LAST component at 0 whatever was written:
# theta2 ~ x on two components was refused as an unknown dpar, and
# theta1 ~ x on three was fitted with an intercept-only theta2.

mref_data <- function() {
  set.seed(109)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  p2 <- stats::plogis(-0.3 + 0.9 * d$x)
  k2 <- stats::rbinom(n, 1, p2) == 1
  d$y <- ifelse(k2, stats::rnorm(n, 2, 0.7), stats::rnorm(n, -1, 0.7))
  d
}

test_that("theta2 ~ x on two components makes component 1 the reference", {
  d <- mref_data()
  f2 <- frm(bf(y ~ 1, theta2 ~ x) + mixture(gaussian(), gaussian()),
            data = d)
  f1 <- frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()),
            data = d)
  fe2 <- fixef(f2)
  fe1 <- fixef(f1)
  expect_true(all(c("theta2_Intercept", "theta2_x") %in% rownames(fe2)))
  expect_false(any(grepl("^theta1", rownames(fe2))))
  # the same model under the other reference: an IDENTITY, the log odds
  # of component 2 against 1 is minus that of 1 against 2, so the
  # coefficients negate and the likelihood is the same one
  expect_equal(unname(fe2[c("theta2_Intercept", "theta2_x"), "Estimate"]),
               -unname(fe1[c("theta1_Intercept", "theta1_x"), "Estimate"]),
               tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f1)),
               tolerance = 1e-8)
  # the response scale of theta2 is component 2's share, which grows with
  # x in these data, and the reference has no predictor to ask about
  p <- as.numeric(frm_linpred(f2, newdata = data.frame(x = c(-1, 0, 1)),
                              dpar = "theta2", type = "response"))
  expect_true(all(diff(p) > 0))
  q <- as.numeric(frm_linpred(f1, newdata = data.frame(x = c(-1, 0, 1)),
                              dpar = "theta1", type = "response"))
  expect_equal(p, 1 - q, tolerance = 1e-5)
  expect_error(frm_linpred(f2, dpar = "theta1"), "theta1",
               class = "frmtmb_error")
})

test_that("three components take a formula for exactly two thetas", {
  d <- mref_data()
  mix3 <- mixture(gaussian(), gaussian(), gaussian())
  expect_error(frm(bf(y ~ 1, theta1 ~ x) + mix3, data = d),
               "all but one mixing proportion", class = "frmtmb_error")
  expect_error(default_prior(bf(y ~ 1, theta1 ~ x) + mix3, data = d),
               "all but one mixing proportion", class = "frmtmb_error")
  expect_error(frm(bf(y ~ 1, theta1 ~ x, theta2 ~ x) +
                     mixture(gaussian(), gaussian()), data = d),
               "all but one mixing proportion", class = "frmtmb_error")
  # theta2 and theta3 written: component 1 is the reference
  tab <- default_prior(bf(y ~ 1, theta2 ~ x, theta3 ~ x) + mix3, data = d)
  expect_setequal(unique(tab$dpar[grepl("^theta", tab$dpar)]),
                  c("theta2", "theta3"))
  # and with no theta formula the mixture is unchanged: brms's simplex
  # theta1 ... thetaK, the last component the reference internally
  f0 <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = d)
  expect_true(all(c("theta1", "theta2") %in% variables(f0)))
})
