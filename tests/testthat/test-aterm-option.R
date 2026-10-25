# frmtmb_register_aterm(option =): an addition term that carries no data
# and sets an argument of the family instead, as brms's bhaz() sets the
# baseline basis of cox() (frmtmb.survival). The tests register names no
# extension uses and restore every registry (helper-registry.R).

# a gaussian whose mean is shifted by a constant the family holds
shifted_family <- function(shift) {
  custom_family(
    "tshifted", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu + shift, dpars$sigma, log = TRUE)
    })
}

option_data <- function() {
  set.seed(41)
  dd <- data.frame(x = stats::rnorm(80))
  dd$y <- stats::rnorm(80, 1 + 0.5 * dd$x)
  dd
}

test_that("an option term gives the fit of the family it returns", {
  local_frmtmb_registries()
  seen <- new.env()
  frmtmb_register_aterm("tshift", option = function(family, call, env) {
    seen$call <- call
    seen$family <- family[["family"]]
    shifted_family(eval(call[[2L]], env))
  })
  dd <- option_data()
  sh <- 2
  a <- frm(bf(y | tshift(sh) ~ x) + shifted_family(0), data = dd)
  b <- frm(bf(y ~ x) + shifted_family(2), data = dd)
  expect_identical(logLik(a), logLik(b))
  expect_identical(a$estimates, b$estimates)
  # the function saw the call as written and the family as given, and
  # read the call's argument in the formula's environment
  expect_identical(seen$call, quote(tshift(sh)))
  expect_identical(seen$family, "tshifted")
  # the term is not data: it is no addition-term value, and new data
  # need not carry it
  spec <- frm(bf(y | tshift(sh) ~ x) + shifted_family(0), data = NULL,
              dry_run = "spec")
  expect_false("tshift" %in% names(spec$responses[[1L]]$aterms))
  expect_identical(frm_linpred(a, newdata = data.frame(x = 0:2)),
                   frm_linpred(b, newdata = data.frame(x = 0:2)))
  # the absent case: without the term, the family is the one written
  expect_false(identical(
    logLik(frm(bf(y ~ x) + shifted_family(0), data = dd)), logLik(a)))
})

test_that("an option term composes with data terms and with mvbf()", {
  # what frm_compat() says of the generic rows: an option term adds
  # nothing to the likelihood, and in mvbf() it sets its own response's
  # family only
  local_frmtmb_registries()
  frmtmb_register_aterm("tshift", option = function(family, call, env) {
    shifted_family(eval(call[[2L]], env))
  })
  dd <- option_data()
  dd$w <- rep(c(1, 2), 40)
  set.seed(42)
  dd$z <- stats::rnorm(80, 0.5 * dd$x)
  a <- frm(bf(y | weights(w) + tshift(2) ~ x) + shifted_family(0),
           data = dd)
  b <- frm(bf(y | weights(w) ~ x) + shifted_family(2), data = dd)
  expect_identical(logLik(a), logLik(b))
  expect_identical(a$estimates, b$estimates)
  m1 <- frm(mvbf(bf(y | tshift(2) ~ x) + shifted_family(0),
                 bf(z ~ x) + shifted_family(0)) + set_rescor(FALSE),
            data = dd)
  m2 <- frm(mvbf(bf(y ~ x) + shifted_family(2),
                 bf(z ~ x) + shifted_family(0)) + set_rescor(FALSE),
            data = dd)
  expect_identical(logLik(m1), logLik(m2))
  rows <- frm_compat("tshift()")
  expect_match(rows$note[rows$feature_b == "weights()"],
               "an option term", fixed = TRUE)
  expect_match(rows$note[rows$feature_b == "mvbf"],
               "sets the family of its own response only", fixed = TRUE)
})

test_that("a family that does not declare the term is refused first", {
  local_frmtmb_registries()
  calls <- 0L
  frmtmb_register_aterm("tshift", option = function(family, call, env) {
    calls <<- calls + 1L
    family
  })
  expect_error(frm(bf(y | tshift(2) ~ x), family = gaussian(), data = NULL,
                   dry_run = "spec"),
               "Addition term `tshift()` is not supported for family",
               fixed = TRUE)
  expect_identical(calls, 0L)
  # a family with no allow-list reaches the function
  frm(bf(y | tshift(2) ~ x) + shifted_family(0), data = NULL,
      dry_run = "spec")
  expect_identical(calls, 1L)
})

test_that("an option term is checked like any other", {
  local_frmtmb_registries()
  frmtmb_register_aterm("tshift", option = function(family, call, env) 1)
  expect_error(frm(bf(y | tshift(2) ~ x) + shifted_family(0), data = NULL,
                   dry_run = "spec"),
               "must return the family the fit uses", fixed = TRUE)
  frmtmb_register_aterm("tshift", option = function(family, call, env) {
    family
  })
  expect_error(frm(bf(y | tshift(1) + tshift(2) ~ x) + shifted_family(0),
                   data = NULL, dry_run = "spec"),
               "Duplicated addition term `tshift()`", fixed = TRUE)
  expect_error(frmtmb_register_aterm("tshift2", option = "shift"),
               "option =) must be NULL or a function", fixed = TRUE)
  # a data term registered beside it is unchanged
  frmtmb_register_aterm("tside")
  spec <- frm(bf(y | tside(b) + tshift(2) ~ x) + shifted_family(0),
              data = NULL, dry_run = "spec")
  expect_identical(names(spec$responses[[1L]]$aterms), "tside")
})
