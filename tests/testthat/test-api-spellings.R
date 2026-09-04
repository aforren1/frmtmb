# The brms argument-naming rule where it has to be dual: a brms-NAMED
# function speaks brms's argument names, but the draws surface SHIPPED
# taking lme4's `re.form`, so both spellings stay live and mean one
# thing. What is asserted here is the compatibility claim itself - the
# alias resolves to the same internal setting, the defaults are the ones
# that were there before the alias existed, and giving both spellings at
# once is refused rather than silently resolved.

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
                       g = factor(rep(1:6, 10)))
      dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                             stats::rnorm(6, 0, 0.5)[dd$g], 1)
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      cache <<- list(dd = dd, fit = fit)
    }
    cache
  }
})

# ---- the surface itself ----------------------------------------------

test_that("the brms-named fit method carries both spellings", {
  # a structural guard, so a method added later cannot quietly ship one
  # spelling: the list is the contract. pp_check() on a FIT is the only
  # brms-named frmtmb method that takes a re_formula; the draws methods
  # that were in this list are asserted the same way, in the same
  # words, by frmtmb.sample.
  fo <- formals(getFromNamespace("pp_check.frmtmb_fit", "frmtmb"))
  expect_true(all(c("re_formula", "re.form") %in% names(fo)))
  ns <- asNamespace("frmtmb")
  expect_true(frmtmb:::is_arg_unset(eval(fo[["re_formula"]], ns)))
  expect_true(frmtmb:::is_arg_unset(eval(fo[["re.form"]], ns)))
})

test_that("giving both spellings is refused on a fit, not resolved", {
  cs <- sp_case()
  skip_if_not_installed("bayesplot")

  expect_error(pp_check(cs$fit, re_formula = NA, re.form = NA),
               "pp_check")
  # the refusal names both spellings and says which one the function is
  # named after, so it can be acted on without reading the manual
  msg <- tryCatch(pp_check(cs$fit, re_formula = NA, re.form = NA),
                  error = conditionMessage)
  expect_match(msg, "re_formula")
  expect_match(msg, "re.form")

  # agreeing values are refused too: the point is that the call did not
  # say which name it meant, not that the two disagreed
  expect_error(pp_check(cs$fit, re_formula = NULL, re.form = NULL),
               "two spellings of ONE setting")
})

test_that("re_formula and re.form give identical pp_check() output", {
  skip_if_not_installed("bayesplot")
  cs <- sp_case()

  set.seed(7)
  p1 <- pp_check(cs$fit, ndraws = 5, re_formula = NA)
  set.seed(7)
  p2 <- pp_check(cs$fit, ndraws = 5, re.form = NA)
  expect_equal(p1$data, p2$data)
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


# ---- defaults, unchanged ---------------------------------------------


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
