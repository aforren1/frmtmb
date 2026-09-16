# The argument-naming rule: brms is the TIEBREAKER. Where
# lme4 or glmmTMB and brms disagree on a name, this package takes brms's
# and DROPS the other one; whether the other one already shipped is not
# a consideration. Two spellings of one setting survive only where brms
# ITSELF carries both, which is `posterior_epred()`,
# `posterior_linpred()`, `posterior_predict()` and `predictive_error()`
# on a `brmsfit`, and those are asserted in frmtmb.sample's
# test-draws-spellings.R. What is asserted here is that the FIT surface
# now speaks brms alone.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

# The fixture no longer samples, and takes no tmbstan skip: the draws
# half of this suite is frmtmb.sample's test-draws-spellings.R, and
# what is left here needs only the fit.
sp_case <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(11)
      dd <- data.frame(x = stats::rnorm(60),
                       g = factor(rep(1:6, 10)), y = 0)
      # this file compares argument SPELLINGS, so the response is drawn
      # from the model rather than restated in rnorm() calls
      # The draw takes its own seed, away from the fixture's: reusing
      # the fixture seed restarts the same random stream that made the
      # covariates, and the residuals come out equal to x.
      dd$y <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                           newparams = list(Intercept = 1, x = 0.5,
                                            sigma = 1,
                                            sd_g__Intercept = 0.5),
                           nsim = 1, seed = 1011)[[1]]
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      cache <<- list(dd = dd, fit = fit)
    }
    cache
  }
})

# ---- the surface itself ----------------------------------------------

test_that("pp_check() on a fit takes re_formula alone", {
  # brms's own pp_check.brmsfit() forwards to prepare_predictions(),
  # whose formals carry `re_formula` and NOT `re.form`, so the alias
  # this method used to accept had no brms precedent behind it.
  fo <- names(formals(getFromNamespace("pp_check.frmtmb_fit", "frmtmb")))
  expect_true("re_formula" %in% fo)
  expect_false("re.form" %in% fo)
})

test_that("the retired lme4 spelling is refused by pp_check(), not passed on", {
  # pp_check() forwards its dots to bayesplot's ppc_* function, which
  # takes dots of its own and would have accepted `re.form` and done
  # nothing with it. The retired spelling is therefore checked here
  # before the forward, which is the only name that is.
  fo <- names(formals(getFromNamespace("pp_check.frmtmb_fit", "frmtmb")))
  expect_false("re.form" %in% fo)
  cs <- sp_case()
  skip_if_not_installed("bayesplot")
  expect_error(pp_check(cs$fit, ndraws = 5, re.form = NA), "re_formula")
  expect_error(predict(cs$fit, re.form = NA), "re.form", fixed = TRUE)
})

test_that("the fit surface speaks brms alone", {
  # the rule itself: predict(), simulate() and frm_bootstrap() take
  # brms's `re_formula` and lme4's `re.form` is gone from all three
  for (nm in c("predict.frmtmb_fit", "simulate.frmtmb_fit")) {
    fo <- names(formals(getFromNamespace(nm, "frmtmb")))
    expect_true("re_formula" %in% fo)
    expect_false("re.form" %in% fo)
  }
  fo <- names(formals(frm_bootstrap))
  expect_true("re_formula" %in% fo)
  expect_false("re.form" %in% fo)
})

# ---- equivalence ------------------------------------------------------


# ---- defaults, unchanged ---------------------------------------------


test_that("pp_check() on a fit still defaults to NA", {
  skip_if_not_installed("bayesplot")
  cs <- sp_case()

  set.seed(10)
  p0 <- pp_check(cs$fit, ndraws = 5)
  set.seed(10)
  pna <- pp_check(cs$fit, ndraws = 5, re_formula = NA)
  expect_equal(p0$data, pna$data)

  # NA simulates new levels, NULL reuses the fitted modes, so the
  # default is identifiable rather than merely asserted
  set.seed(10)
  pnull <- pp_check(cs$fit, ndraws = 5, re_formula = NULL)
  expect_false(isTRUE(all.equal(p0$data, pnull$data)))
})

# ---- both spellings at once ------------------------------------------


# ---- set_rescor(), the other name divergence the diff found ----------

test_that("set_rescor() takes brms's rescor and keeps rescor_value", {
  expect_identical(set_rescor(TRUE), set_rescor(rescor = TRUE))
  expect_identical(set_rescor(rescor = FALSE),
                   set_rescor(rescor_value = FALSE))
  expect_true(set_rescor()$rescor)          # the default is still TRUE
  expect_false(set_rescor(FALSE)$rescor)
  expect_error(set_rescor(rescor = TRUE, rescor_value = TRUE),
               "two spellings of ONE setting")
  expect_error(set_rescor(rescor = "yes"), "must be TRUE or FALSE")
})

# --- one coefficient, one name ----------------------------------------
# fixef() is a LIST keyed by dpar whose entries are named by design
# COLUMN, so unlist() composes base R's own "mu.x" - a name no surface
# of this package uses. vcov(), confint(), hypothesis(), par_template()
# and `start` all spell the same coefficient "sigma_x", with the
# LOCATION parameter's own coefficients left unprefixed. flatten = TRUE
# is that spelling; the unlist() composite is pinned below so it cannot
# drift into looking canonical.

test_that("fixef(flatten = TRUE) names coefficients as vcov() names its rows", {
  cs <- sp_case()
  f <- cs$fit
  expect_setequal(names(fixef(f, flatten = TRUE)), rownames(vcov(f)))
  # the location parameter is NOT prefixed, here or there
  expect_true("x" %in% names(fixef(f, flatten = TRUE)))
  expect_false("mu.x" %in% names(fixef(f, flatten = TRUE)))
  # the values are the same coefficients, just differently keyed
  expect_equal(unname(fixef(f, flatten = TRUE)[["x"]]),
               unname(fixef(f)$mu[["x"]]))
  # every flattened name addresses a confint() row
  expect_true(all(names(fixef(f, flatten = TRUE)) %in% rownames(confint(f))))
})

test_that("a second dpar keeps the two spellings apart", {
  set.seed(4)
  dd <- data.frame(x = stats::rnorm(120))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x, exp(0.2 + 0.1 * dd$x))
  f <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
  expect_setequal(names(fixef(f, flatten = TRUE)), rownames(vcov(f)))
  expect_true("sigma_x" %in% names(fixef(f, flatten = TRUE)))
  # the default shape is untouched: a list keyed by dpar, bare columns
  expect_named(fixef(f), c("mu", "sigma"))
  expect_named(fixef(f)$sigma, c("(Intercept)", "x"))
})

test_that("unlist(fixef()) stays base R's composite, deliberately", {
  # NOT the canonical spelling, and not made into one: 19 files in this
  # monorepo index literal "dpar.column" names out of unlist(fixef()).
  # Pinned so a later change has to be a decision.
  cs <- sp_case()
  expect_true("mu.x" %in% names(unlist(fixef(cs$fit))))
  expect_false("mu.x" %in% rownames(vcov(cs$fit)))
})

test_that("the bootstrap statistic is named like the covariance it is compared with", {
  cs <- sp_case()
  bs <- frm_bootstrap(cs$fit, nsim = 3, seed = 1)
  expect_setequal(names(bs$t0), rownames(vcov(cs$fit)))
  expect_equal(unname(bs$t0), unname(fixef(cs$fit, flatten = TRUE)))
})

test_that("flatten is a flag", {
  cs <- sp_case()
  expect_error(fixef(cs$fit, flatten = "yes"), "flatten")
})
