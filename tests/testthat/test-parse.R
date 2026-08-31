test_that("bf() builds a formula object and + attaches a family", {
  f <- bf(y ~ x + (x | g))
  expect_s3_class(f, "frmtmb_formula")
  expect_null(f$family)
  f2 <- f + gaussian()
  expect_s3_class(f2$family, "frmtmb_family")
  expect_identical(f2$family$family, "gaussian")
})

test_that("bf() rejects unsupported features with clear errors", {
  expect_error(bf(y ~ a * exp(b * x), nl = TRUE), "parameter formula")
  expect_error(frm(bf(y ~ x), data.frame(y = 1:3, x = 1:3)),
               "No family")
})

test_that("bf() collects dpar formulas and constants", {
  f <- bf(y ~ x, sigma ~ z + (1 | g))
  expect_named(f$pforms, "sigma")
  f2 <- bf(y ~ x, sigma = 1.5)
  expect_identical(f2$pfix$sigma, 1.5)
  spec <- frm(bf(y ~ x, sigma ~ z + (1 | g)) + gaussian(),
                 data = NULL, dry_run = "spec")
  dp <- spec$responses[[1]]$dpars$sigma
  expect_identical(deparse1(dp$fixed), "~z")
  expect_length(dp$re, 1)
})

test_that("parse_spec builds the IR for random-effect terms", {
  spec <- frm(bf(y ~ x + (x | g)) + poisson(),
                 data = NULL, dry_run = "spec")
  expect_s3_class(spec, "frmtmb_spec")
  r <- spec$responses[[1]]
  expect_identical(r$resp_name, "y")
  expect_identical(r$family$family, "poisson")
  expect_length(r$dpars$mu$re, 1)
  expect_identical(r$dpars$mu$re[[1]]$covstruct, "us")
  expect_identical(deparse1(r$dpars$mu$fixed), "~x")
})

test_that("double-bar and diag() terms parse to independent structures", {
  spec <- frm(bf(y ~ x + (x || g)) + poisson(),
                 data = NULL, dry_run = "spec")
  re <- spec$responses[[1]]$dpars$mu$re
  expect_length(re, 2)

  spec2 <- frm(bf(y ~ x + diag(x | g)) + poisson(),
                  data = NULL, dry_run = "spec")
  re2 <- spec2$responses[[1]]$dpars$mu$re
  expect_identical(re2[[1]]$covstruct, "diag")
})

test_that("aterms parse and unsupported aterms error", {
  spec <- frm(bf(y | trials(n) ~ x) + binomial(),
                 data = NULL, dry_run = "spec")
  expect_named(spec$responses[[1]]$aterms, "trials")
  expect_error(bf_spec <- frm(bf(y | se(s2) ~ x) + gaussian(),
                                 data = NULL, dry_run = "spec"),
               "se")
})

test_that("gaussian gets an intercept-only sigma dpar", {
  spec <- frm(bf(y ~ x) + gaussian(), data = NULL, dry_run = "spec")
  dp <- spec$responses[[1]]$dpars
  expect_named(dp, c("mu", "sigma"))
  expect_identical(deparse1(dp$sigma$fixed), "~1")
  expect_identical(dp$sigma$link$name, "log")
})

test_that("unsupported covariance structures error at parse time", {
  expect_error(frm(bf(y ~ toep(x | g)) + gaussian(),
                      data = NULL, dry_run = "spec"),
               "not supported yet")
})
