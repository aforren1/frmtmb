# brms names a mi() predictor's coefficient bsp_<resp>_mi<x>, in the
# special-term vector it shares with mo(), and lists it after every b_
# coefficient (brms 2.23.0, dev/correct-log/brms-mi.txt, script
# dev/correct-brms-probe.R, data seed 108):
#   variables: b_y_Intercept b_xm_Intercept b_y_z b_xm_z bsp_y_mixm ...
#   fixef rows: y_Intercept xm_Intercept y_z xm_z y_mixm
#   interaction: bsp_y_mixm bsp_y_mixm:z; fixef ... y_mixm y_mixm:z
#   prior rows: class b, coef mixm, resp y (unchanged by the name)
#   hypothesis(fit, "y_mixm > 0"): refused, 'b_y_mixm' not found
# Through 0.61.0 frmtmb named it b_y_mixm and ordered the fixef rows by
# design column.

mib_data <- function() {
  set.seed(108)
  n <- 80
  d <- data.frame(z = stats::rnorm(n))
  d$xm <- 0.5 * d$z + stats::rnorm(n)
  d$y <- 1 + 0.7 * d$xm - 0.3 * d$z + stats::rnorm(n)
  d$xm[sample.int(n, 12)] <- NA
  d
}

mib_fit <- function() {
  frm(bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z) + set_rescor(FALSE) +
        gaussian(), data = mib_data())
}

test_that("a mi() coefficient is bsp_, as in brms", {
  fit <- mib_fit()
  v <- variables(fit)
  expect_true("bsp_y_mixm" %in% v)
  expect_false("b_y_mixm" %in% v)
  expect_identical(rownames(fixef(fit)),
                   c("y_Intercept", "xm_Intercept", "y_z", "xm_z", "y_mixm"))
  expect_identical(rownames(summary(fit)$fixed), rownames(fixef(fit)))
  expect_identical(rownames(vcov(fit)), rownames(fixef(fit)))
})

# Its own block: on 0.61.0 the refusal below is not a frmtmb_error but
# "Some parameters cannot be found in the model", so testthat ended the
# block there and what follows never ran on the base arm.
test_that("hypothesis() takes the bsp_ name and refuses the b_ one", {
  fit <- mib_fit()
  # brms's default class = "b" puts b_ in front, which no longer names
  # this coefficient, so brms refuses the bare name and so does frmtmb
  expect_error(hypothesis(fit, "y_mixm > 0"), "b_y_mixm",
               class = "frmtmb_error")
  h <- hypothesis(fit, "bsp_y_mixm > 0", class = NULL)
  expect_equal(h$hypothesis$Estimate,
               unname(fixef(fit)["y_mixm", "Estimate"]))
})

test_that("the mi() prior row keeps brms's class and coef", {
  d <- mib_data()
  # the prior row is brms's: class b, the coef without the prefix
  tab <- default_prior(bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z) +
                         set_rescor(FALSE), d)
  expect_true(any(tab$class == "b" & tab$coef == "mixm" & tab$resp == "y"))
})

test_that("a mi() interaction lists its main effect first, as in brms", {
  d <- mib_data()
  fit <- frm(bf(y ~ mi(xm) * z) + bf(xm | mi() ~ z) + set_rescor(FALSE) +
               gaussian(), data = d)
  v <- variables(fit)
  expect_true(all(c("bsp_y_mixm", "bsp_y_mixm:z") %in% v))
  # brms enumerates its special terms in terms() order, main effects
  # first, and variables() follows it here because the design columns do
  expect_lt(match("bsp_y_mixm", v), match("bsp_y_mixm:z", v))
  expect_identical(rownames(fixef(fit)),
                   c("y_Intercept", "xm_Intercept", "y_z", "xm_z", "y_mixm",
                     "y_mixm:z"))
})
