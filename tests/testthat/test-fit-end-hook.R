# The fit-end hook: a family's own check, run once when a fit finishes,
# and what happens when it fails.
#
# `logLik()` reads the optimizer's own value, so a family whose
# likelihood is floored or degenerate somewhere had no way to say so.
# `post$fit_check` is that way. The property this file exists for is the
# second one: a check runs AFTER frm() has done all the work, so a check
# that throws must not take the fit with it.

hook_family <- function(body) {
  frmtmb_family(
    "hooked", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) stats::sd(y)),
    post = list(mean_fn = function(dpars, aterms) dpars[["mu"]],
                fit_check = body))
}

hook_data <- function(seed = 8081) {
  set.seed(seed)
  d <- data.frame(x = rnorm(80))
  d$y <- rnorm(80, 1 + 2 * d$x, 0.5)
  d
}

test_that("a family's fit_check hook runs when the fit finishes", {
  skip_on_cran()
  d <- hook_data()
  seen <- new.env(parent = emptyenv())
  fam <- hook_family(function(fit, resp) {
    seen$resp <- resp
    seen$fitted <- inherits(fit, "frmtmb_fit") && !is.null(fit$opt)
    invisible(NULL)
  })
  fit <- frm(bf(y ~ x), d, fam)
  expect_identical(seen$resp, "y")
  # it sees a FINISHED fit, which is the whole point: nothing else can
  expect_true(seen$fitted)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("a hook that warns reaches the user and keeps the fit", {
  skip_on_cran()
  d <- hook_data()
  fam <- hook_family(function(fit, resp) {
    warning("hooked: this family is unhappy about ", resp, call. = FALSE)
  })
  expect_warning(fit <- frm(bf(y ~ x), d, fam),
                 "this family is unhappy about y")
  expect_s3_class(fit, "frmtmb_fit")
  expect_equal(unname(fixef(fit)$mu[["x"]]), 2, tolerance = 0.1)
})

test_that("a hook that THROWS does not destroy the fit", {
  skip_on_cran()
  d <- hook_data()
  fam <- hook_family(function(fit, resp) {
    stop("hooked: deliberate failure inside the check", call. = FALSE)
  })
  # the fit is returned, the failure is reported, and the message names
  # the family and the hook so the reader knows whose check broke
  expect_warning(fit <- frm(bf(y ~ x), d, fam), "post\\$fit_check hook failed")
  w <- tryCatch(frm(bf(y ~ x), d, fam),
                warning = function(e) conditionMessage(e))
  expect_true(grepl("hooked", w, fixed = TRUE))
  expect_true(grepl("deliberate failure inside the check", w, fixed = TRUE))
  expect_true(grepl("complete and unaffected", w, fixed = TRUE))

  expect_s3_class(fit, "frmtmb_fit")
  expect_equal(unname(fixef(fit)$mu[["x"]]), 2, tolerance = 0.1)
  expect_true(is.finite(as.numeric(logLik(fit))))
  # and the fit is the SAME fit a family with no hook would have given
  plain <- frm(bf(y ~ x), d, hook_family(NULL))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(plain)),
               tolerance = 1e-10)
})

test_that("the guard covers core's own coverage report too", {
  skip_on_cran()
  # `ps_coverage_warning()` evaluates a nonlinear body to find its
  # arguments, which is arbitrary code as much as a family's hook is, so
  # it sits inside the same guard.
  set.seed(8082)
  n_id <- 20
  d <- expand.grid(t = seq(0, 1, length.out = 8), id = factor(1:n_id))
  d$y <- 2 + sin(2 * pi * d$t) + rnorm(nrow(d), 0, 0.2)
  fit <- frm(bf(y ~ lev + ps(t, k = 8, pad = 0.3), lev ~ 1 + (1 | id),
                nl = TRUE), d, gaussian())
  expect_s3_class(fit, "frmtmb_fit")

  # the guard, exercised directly on a fit whose family carries a
  # throwing hook: the fit object is real, only the hook is planted
  bad <- fit
  bad$spec$responses[[1]]$family$post$fit_check <- function(f, r) {
    stop("boom", call. = FALSE)
  }
  bad$spec$responses[[1]]$family$family <- "planted"
  expect_warning(frmtmb:::fit_end_checks(bad), "post\\$fit_check hook failed")
  expect_warning(frmtmb:::fit_end_checks(bad), "planted")
  # and it returns rather than propagating
  expect_null(suppressWarnings(frmtmb:::fit_end_checks(bad)))
})
