# The brms argument-naming rule where it has to be dual: a brms-NAMED
# function speaks brms's argument names, but the draws surface SHIPPED
# taking lme4's `re.form`, so both spellings stay live and mean one
# thing. What is asserted here is the compatibility claim itself - the
# alias resolves to the same internal setting, the defaults are the ones
# that were there before the alias existed, and giving both spellings at
# once is refused rather than silently resolved.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

skip_sampler <- function() {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
}

sp_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(11)
      dd <- data.frame(x = stats::rnorm(60),
                       g = factor(rep(1:6, 10)))
      dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                             stats::rnorm(6, 0, 0.5)[dd$g], 1)
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 3)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

# ---- the surface itself ----------------------------------------------

test_that("every brms-named draws method carries both spellings", {
  # a structural guard, so a method added later to this family cannot
  # quietly ship one spelling: the list is the contract
  dual <- c("posterior_epred.frmtmb_draws", "posterior_linpred.frmtmb_draws",
            "posterior_predict.frmtmb_draws",
            "predictive_interval.frmtmb_draws",
            "predictive_error.frmtmb_draws",
            "pp_check.frmtmb_draws", "pp_check.frmtmb_fit")
  for (nm in dual) {
    fo <- formals(getFromNamespace(nm, "frmtmb"))
    expect_true(all(c("re_formula", "re.form") %in% names(fo)),
                label = paste0(nm, " has both spellings"))
    # both default to the "not supplied" marker, so either alone is a
    # setting and neither NULL nor NA is mistaken for one
    ns <- asNamespace("frmtmb")
    expect_true(frmtmb:::is_arg_unset(eval(fo[["re_formula"]], ns)))
    expect_true(frmtmb:::is_arg_unset(eval(fo[["re.form"]], ns)))
  }
})

test_that("the fit surface keeps lme4's spelling alone", {
  # the other half of the rule: predict(), simulate() and frm_bootstrap()
  # are frmtmb's own, not brms's, so they do NOT gain re_formula
  for (nm in c("predict.frmtmb_fit", "simulate.frmtmb_fit")) {
    fo <- names(formals(getFromNamespace(nm, "frmtmb")))
    expect_true("re.form" %in% fo)
    expect_false("re_formula" %in% fo)
  }
  fo <- names(formals(frm_bootstrap))
  expect_true("re.form" %in% fo)
  expect_false("re_formula" %in% fo)
})

# ---- equivalence ------------------------------------------------------

test_that("re_formula and re.form give identical draws-method output", {
  cs <- sp_case()
  ds <- cs$ds
  nd <- data.frame(x = c(-1, 0, 1),
                   g = factor(1, levels = levels(cs$dd$g)))

  for (rf in list(NA, NULL, ~0)) {
    expect_equal(posterior_epred(ds, newdata = nd, re_formula = rf,
                                 ndraws = 12),
                 posterior_epred(ds, newdata = nd, re.form = rf,
                                 ndraws = 12))
    expect_equal(posterior_linpred(ds, newdata = nd, re_formula = rf,
                                   ndraws = 12),
                 posterior_linpred(ds, newdata = nd, re.form = rf,
                                   ndraws = 12))
  }

  # the simulating methods need the same seed to be comparable at all
  set.seed(4)
  a <- posterior_predict(ds, newdata = nd, re_formula = NA, ndraws = 12)
  set.seed(4)
  b <- posterior_predict(ds, newdata = nd, re.form = NA, ndraws = 12)
  expect_equal(a, b)

  set.seed(5)
  a <- predictive_interval(ds, re_formula = NA, ndraws = 12)
  set.seed(5)
  b <- predictive_interval(ds, re.form = NA, ndraws = 12)
  expect_equal(a, b)

  set.seed(6)
  a <- predictive_error(ds, re_formula = NA, ndraws = 12)
  set.seed(6)
  b <- predictive_error(ds, re.form = NA, ndraws = 12)
  expect_equal(a, b)
})

test_that("re_formula and re.form give identical pp_check() output", {
  skip_if_not_installed("bayesplot")
  cs <- sp_case()

  set.seed(7)
  p1 <- pp_check(cs$fit, ndraws = 5, re_formula = NA)
  set.seed(7)
  p2 <- pp_check(cs$fit, ndraws = 5, re.form = NA)
  expect_equal(p1$data, p2$data)

  set.seed(8)
  d1 <- pp_check(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(8)
  d2 <- pp_check(cs$ds, ndraws = 5, re.form = NA)
  expect_equal(d1$data, d2$data)

  # non-vacuity: the switch takes effect IN-SAMPLE (review finding: it
  # used to be consulted only under newdata), so NA must differ from
  # the default NULL wherever the model has random effects, under the
  # same RNG seed
  set.seed(9)
  dNA <- posterior_predict(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(9)
  dNU <- posterior_predict(cs$ds, ndraws = 5)
  expect_false(isTRUE(all.equal(dNA, dNU)))
  set.seed(9)
  dNA2 <- posterior_predict(cs$ds, ndraws = 5, re.form = NA)
  expect_equal(dNA, dNA2)
  set.seed(9)
  eNA <- predictive_error(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(9)
  eNU <- predictive_error(cs$ds, ndraws = 5)
  expect_false(isTRUE(all.equal(eNA, eNU)))
  set.seed(9)
  pNA <- pp_check(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(9)
  pNU <- pp_check(cs$ds, ndraws = 5)
  expect_false(isTRUE(all.equal(pNA$data, pNU$data)))
})

# ---- defaults, unchanged ---------------------------------------------

test_that("the draws methods still default to NULL", {
  cs <- sp_case()
  ds <- cs$ds

  # the pinned pre-change default of every draws method: condition on
  # each draw's own random effects
  expect_equal(posterior_epred(ds, ndraws = 12),
               posterior_epred(ds, re.form = NULL, ndraws = 12))
  expect_equal(posterior_epred(ds, ndraws = 12),
               posterior_epred(ds, re_formula = NULL, ndraws = 12))
  expect_equal(posterior_linpred(ds, ndraws = 12),
               posterior_linpred(ds, re.form = NULL, ndraws = 12))

  # and NULL is a different quantity from NA, so the assertion above is
  # not vacuous: a default silently flipped to NA would show up here
  expect_false(isTRUE(all.equal(
    posterior_epred(ds, ndraws = 12),
    posterior_epred(ds, re_formula = NA, ndraws = 12))))

  set.seed(9)
  a <- predictive_error(ds, ndraws = 12)
  set.seed(9)
  b <- predictive_error(ds, re.form = NULL, ndraws = 12)
  expect_equal(a, b)
})

test_that("pp_check() on a fit still defaults to NA", {
  skip_if_not_installed("bayesplot")
  cs <- sp_case()

  set.seed(10)
  p0 <- pp_check(cs$fit, ndraws = 5)
  set.seed(10)
  pna <- pp_check(cs$fit, ndraws = 5, re.form = NA)
  expect_equal(p0$data, pna$data)

  # NA simulates new levels, NULL reuses the fitted modes, so the
  # default is identifiable rather than merely asserted
  set.seed(10)
  pnull <- pp_check(cs$fit, ndraws = 5, re.form = NULL)
  expect_false(isTRUE(all.equal(p0$data, pnull$data)))
})

# ---- both spellings at once ------------------------------------------

test_that("giving both spellings is refused, not resolved", {
  cs <- sp_case()
  ds <- cs$ds

  expect_error(posterior_epred(ds, re_formula = NA, re.form = NA),
               "two spellings of ONE setting")
  expect_error(posterior_linpred(ds, re_formula = NA, re.form = NA),
               "posterior_linpred\\(\\)")
  expect_error(posterior_predict(ds, re_formula = NA, re.form = NA),
               "posterior_predict\\(\\)")
  expect_error(predictive_interval(ds, re_formula = NA, re.form = NA),
               "predictive_interval\\(\\)")
  expect_error(predictive_error(ds, re_formula = NA, re.form = NA),
               "predictive_error\\(\\)")
  expect_error(pp_check(ds, re_formula = NA, re.form = NA),
               "pp_check\\(\\)")
  expect_error(pp_check(cs$fit, re_formula = NA, re.form = NA),
               "pp_check\\(\\)")

  # the refusal names both spellings and says which one the function is
  # named after, so it can be acted on without reading the manual
  msg <- tryCatch(posterior_epred(ds, re_formula = NA, re.form = NA),
                  error = conditionMessage)
  expect_match(msg, "`re_formula`")
  expect_match(msg, "`re.form`")

  # agreeing values are refused too: the point is that the call did not
  # say which name it meant, not that the two disagreed
  expect_error(posterior_epred(ds, re_formula = NULL, re.form = NULL),
               "two spellings of ONE setting")
})

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
