# The sampler half of lane wt-adefects: `log_lik()` is frmtmb's generic
# now, re-exported here, so that a user who has not sampled is told to
# sample rather than told the function does not exist
# (dev/adefects-findings.md, D7). This file was SEEN FAILING against the
# base build, where frmtmb.sample defined its own generic and frmtmb had
# none.

test_that("log_lik is frmtmb's generic, re-exported, not a rival one", {
  expect_identical(frmtmb.sample::log_lik, frmtmb::log_lik)
  # the DIRECTIVE, not the resolved value: with rstantools loaded both
  # packages' bindings resolve to rstantools' generic, so comparing the
  # values alone would hold on a build that still defines a rival here.
  # Read with parseNamespaceFile(), never with a grep: roxygen writes
  # one multi-name importFrom() block per package, which a
  # one-name-per-line grep cannot see (R/generic-owners.R records the
  # same trap for nlme).
  ns <- parseNamespaceFile("frmtmb.sample",
                           dirname(find.package("frmtmb.sample")))
  from_frmtmb <- unlist(lapply(ns$imports, function(x) {
    if (is.list(x) && identical(x[[1L]], "frmtmb")) x[[2L]] else character(0)
  }))
  expect_true("log_lik" %in% from_frmtmb)

  # it is NOT in this package's own owner table any more: a re-export
  # inherits frmtmb's active binding rather than installing a second one
  own <- get("sample_generic_owners", envir = asNamespace("frmtmb.sample"))
  expect_false("log_lik" %in% names(own))

  # the draws method is still registered on the owner's generic too, so
  # it is reached once rstantools takes the binding
  m <- parseNamespaceFile("frmtmb.sample",
                          dirname(find.package("frmtmb.sample")))$S3methods
  row <- m[, 1] == "log_lik" & m[, 2] == "frmtmb_draws"
  expect_equal(sum(row), 2L)
  expect_true(any(is.na(m[row, 4])))
  expect_true(any(m[row, 4] %in% "rstantools"))
})

test_that("log_lik on a fit refuses with this package attached", {
  # the defect this closes: with frmtmb.sample attached, `log_lik(fit)`
  # reached a generic whose method table had no frmtmb_fit entry
  set.seed(20260917)
  dd <- data.frame(x = rnorm(40))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(y ~ x, dd)
  expect_error(log_lik(fit), class = "frmtmb_error")
  expect_error(log_lik(fit), "Sample first")
})

# Punch round 1, MINOR 3. Both methods used to swallow `allow_new_levels`
# in `...`, and the caller then met frmtmb's own refusal, whose remedy
# named that same argument. brms FORWARDS it: measured on
# brmsfit_example1 with `visit` removed from newdata,
# posterior_epred() and posterior_predict() both answer 25 x 3 under
# allow_new_levels = TRUE and both die on base R's "object 'visit' not
# found" without it (dev/adefects-findings.md, punch round 1).
ad_draws <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(20260918)
      dd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
      dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                             stats::rnorm(6, 0, 0.5)[dd$g], 1)
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      cache <<- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 3)))
    }
    cache
  }
})

test_that("the draws methods forward allow_new_levels, as brms does", {
  skip_on_cran()
  ds <- ad_draws()
  nd <- data.frame(x = c(-1, 0, 1))

  for (f in list(posterior_epred, posterior_predict)) {
    # the grouping column may be left out entirely
    ep <- f(ds, newdata = nd, allow_new_levels = TRUE, ndraws = 3)
    expect_equal(dim(ep), c(3L, 3L))
    # and without it the call is refused by the D6 message, which names
    # the argument the method now HAS
    expect_error(f(ds, newdata = nd, ndraws = 3), class = "frmtmb_error")
    expect_error(f(ds, newdata = nd, ndraws = 3), "allow_new_levels = TRUE")
    # every OTHER name in the dots is refused rather than ignored
    expect_error(f(ds, newdata = nd, ndraws = 3, not_an_argument = 2),
                 "not_an_argument")
  }
})
