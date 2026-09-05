# frm(dry_run = "objective") and the frmtmb_unfitted object it returns.
#
# These were written beside the sampling tests, because assembling an
# objective without optimizing it is what the formula route of
# frm_sample() needs. The feature is frm()s, so the tests stayed when
# the sampling surface left.

# The fixture came from test-sample-direct.R with these blocks. Nothing
# here checks a likelihood, so the response comes from the model that is
# about to be assembled: the fixture states that model once instead of
# restating it in rnorm() calls that can drift from the formula below.
sd_data <- function(seed = 9, n = 60L, ng = 6L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)), y = 0)
  # The draw takes its own seed, away from the fixture's: reusing
  # the fixture seed restarts the same random stream that made the
  # covariates, and the residuals come out equal to x.
  dd$y <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                       newparams = list(Intercept = 1, x = 0.5, sigma = 1,
                                        sd_g__Intercept = 0.5),
                       nsim = 1, seed = seed + 1000L)[[1]]
  dd
}

test_that("dry_run = 'objective' stops before the optimizer", {
  dd <- sd_data()
  uf <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd, dry_run = "objective")
  expect_s3_class(uf, "frmtmb_unfitted")
  expect_s3_class(uf, "frmtmb_fit")
  expect_null(uf$opt)
  # the tape is real and evaluable, which is all the sampler needs
  expect_true(is.finite(uf$obj$fn(uf$obj$par)))
  # and the frame is a normal frame, so the draws surface has its
  # structure
  expect_equal(stats::nobs(uf), nrow(dd))
  expect_length(uf$frame$re_blocks, 1L)

  # quadrature has no unfitted form: its tape is calibrated at an optimum
  dp <- dd
  dp$y <- stats::rpois(nrow(dp), exp(0.3 + 0.4 * dp$x))
  expect_error(frm(bf(y ~ x + (1 | g)) + poisson(), dp,
                   quadrature = TRUE, dry_run = "objective"),
               "no unfitted form")
})

test_that("methods needing an ML quantity refuse on an unfitted object", {
  dd <- sd_data()
  uf <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd, dry_run = "objective")
  for (f in list(function() summary(uf), function() stats::vcov(uf),
                 function() stats::confint(uf), function() stats::logLik(uf),
                 function() stats::AIC(uf), function() fixef(uf),
                 function() ranef(uf), function() VarCorr(uf),
                 function() stats::predict(uf), function() stats::fitted(uf),
                 function() stats::residuals(uf),
                 function() stats::simulate(uf), function() print(uf))) {
    expect_error(f(), "needs a fitted model")
  }
  # what does NOT need one keeps working: the model description
  expect_equal(stats::nobs(uf), nrow(dd))
  expect_s3_class(stats::formula(uf), "formula")
  expect_equal(stats::family(uf)$family, "gaussian")
})
