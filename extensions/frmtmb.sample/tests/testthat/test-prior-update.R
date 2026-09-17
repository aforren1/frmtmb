# A prior on the sampling call against the prior a MAP fit carries: the
# call's specification replaces the stored one for its slot, as brms's
# update(prior = ) replaces a row, and other slots stack.

test_that("a call's specification replaces the fit's for the same slot", {
  stored <- frmtmb::set_prior("normal(0, 1)", class = "b", coef = "x") +
    frmtmb::set_prior("exponential(1)", class = "sd")
  call <- frmtmb::set_prior("", class = "b", coef = "x", lb = 0)
  st <- frmtmb.sample:::prior_stack(stored, call)
  specs <- unclass(st$pl)
  on_x <- Filter(function(s) identical(s$coef, "x"), specs)
  # one specification for b_x, the call's, and its bound without the
  # stored density, which brms's replaced row would not keep either
  expect_length(on_x, 1L)
  expect_null(on_x[[1L]]$dist)
  expect_identical(on_x[[1L]]$lb, 0)
  # the stored sd prior, another slot, is kept
  expect_true(any(vapply(specs, function(s) identical(s$class, "sd"),
                         TRUE)))
  # and the stacked prior is not refused as a duplicate
  expect_silent(frmtmb::check_prior_slots(st$pl))
  # a density on the call replaces a density the same way
  st2 <- frmtmb.sample:::prior_stack(
    stored, frmtmb::set_prior("normal(0, 9)", class = "b", coef = "x"))
  on_x2 <- Filter(function(s) identical(s$coef, "x"), unclass(st2$pl))
  expect_length(on_x2, 1L)
  expect_identical(on_x2[[1L]]$dist$scale, 9)
})
