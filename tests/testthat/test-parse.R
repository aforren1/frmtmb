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
  # A missing family is no longer an error: it defaults to gaussian,
  # the brms / lme4 / glmmTMB convention. The internal guard remains
  # for a spec built without going through as_bform().
  expect_error(parse_spec(bf(y ~ x)), "No family")
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
  # se() parses since v0.14 and maps sigma out
  spec_se <- frm(bf(y | se(s2) ~ x) + gaussian(),
                 data = NULL, dry_run = "spec")
  expect_true("se" %in% names(spec_se$responses[[1]]$aterms))
  expect_equal(spec_se$responses[[1]]$dpars$sigma$constant, 1)
  expect_error(frm(bf(y | dec(d) ~ x) + gaussian(),
                   data = NULL, dry_run = "spec"),
               "dec")
})

test_that("gaussian gets an intercept-only sigma dpar", {
  spec <- frm(bf(y ~ x) + gaussian(), data = NULL, dry_run = "spec")
  dp <- spec$responses[[1]]$dpars
  expect_named(dp, c("mu", "sigma"))
  expect_identical(deparse1(dp$sigma$fixed), "~1")
  expect_identical(dp$sigma$link$name, "log")
})

test_that("unsupported covariance structures error at parse time", {
  expect_error(frm(bf(y ~ propto(x | g)) + gaussian(),
                      data = NULL, dry_run = "spec"),
               "not supported yet")
  # rr parses since v0.16 (default rank 2)
  spec <- frm(bf(y ~ rr(x + 0 | g, d = 2)) + gaussian(),
              data = NULL, dry_run = "spec")
  expect_identical(spec$responses[[1]]$dpars$mu$re[[1]]$covstruct, "rr")
  expect_identical(spec$responses[[1]]$dpars$mu$re[[1]]$rank, 2L)
})

# --- the addition-term registry (frmtmb.ddm finding 1) ----------------

test_that("the refusal names vint()/vreal() and the registration seam", {
  # an extension author reading this message has to learn that the set
  # is not closed and what to reach for meanwhile
  err <- tryCatch(frm(bf(y | dec(bound) ~ x) + gaussian(), data = NULL,
                      dry_run = "spec"),
                  error = conditionMessage)
  expect_match(err, "vint() for integers", fixed = TRUE)
  expect_match(err, "frmtmb_register_aterm()", fixed = TRUE)
})

test_that("a registered aterm parses, coerces and reaches the density", {
  on.exit(frmtmb_aterm_registry$reg[["dec"]] <- NULL, add = TRUE)
  frmtmb_register_aterm("dec", arity = 1,
                        coerce = function(x) as.integer(factor(x)) - 1L)

  spec <- frm(bf(y | dec(bound) ~ x) + gaussian(), data = NULL,
              dry_run = "spec")
  expect_identical(deparse1(spec$responses[[1]]$aterms$dec), "bound")

  # not `d`: an aterm variable that newdata omits falls back to the
  # formula environment, so a one-letter name is a test that passes or
  # fails on what some other file left in scope
  set.seed(31)
  dd <- data.frame(x = stats::rnorm(120),
                   bound = factor(rep(c("lower", "upper"), 60)))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x +
                         0.8 * (dd$bound == "upper"), 1)

  seen <- NULL
  fam <- custom_family(
    "dec_reader", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu + 0.8 * aterms$dec, dpars$sigma,
                  log = TRUE)
    },
    required_aterms = "dec",
    post = list(mean_fn = function(dpars, aterms) {
      dpars$mu + 0.8 * aterms$dec
    }))
  fit <- frm(bf(y | dec(bound) ~ x) + fam, data = dd)
  # the factor arrived as the 0/1 the coercion produced, not as level
  # codes and not as a factor
  seen <- fit$frame$aterm_values$y$dec
  expect_identical(sort(unique(seen)), c(0, 1))
  expect_identical(seen[1:2], c(0, 1))
  # and the same coercion runs on newdata
  nd <- data.frame(x = c(0, 0),
                   bound = factor(c("upper", "lower"),
                                  levels = levels(dd$bound)))
  pr <- predict(fit, newdata = nd, type = "response")
  expect_gt(pr[1], pr[2])
})

test_that("a registered aterm is required on newdata that omits it", {
  on.exit(frmtmb_aterm_registry$reg[["dec"]] <- NULL, add = TRUE)
  frmtmb_register_aterm("dec")
  set.seed(32)
  dd <- data.frame(x = stats::rnorm(60), bound = rep(0:1, 30))
  dd$y <- stats::rnorm(60, dd$x + dd$bound)
  fam <- custom_family(
    "dec_reader2", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu + aterms$dec, dpars$sigma, log = TRUE)
    },
    post = list(mean_fn = function(dpars, aterms) dpars$mu + aterms$dec))
  fit <- frm(bf(y | dec(bound) ~ x) + fam, data = dd)
  expect_error(predict(fit, newdata = data.frame(x = 0),
                       type = "response"),
               "could not be evaluated on newdata")
})

test_that("a registered aterm of arity above one numbers its arguments", {
  on.exit(frmtmb_aterm_registry$reg[["pair"]] <- NULL, add = TRUE)
  frmtmb_register_aterm("pair", arity = 2)
  spec <- frm(bf(y | pair(a, b) ~ x) + gaussian(), data = NULL,
              dry_run = "spec")
  expect_named(spec$responses[[1]]$aterms, c("pair1", "pair2"))
  expect_error(frm(bf(y | pair(a) ~ x) + gaussian(), data = NULL,
                   dry_run = "spec"),
               "takes 2 arguments, not 1")
  expect_error(frm(bf(y | pair(a, b) + pair(c, d) ~ x) + gaussian(),
                   data = NULL, dry_run = "spec"),
               "Duplicated addition term")
})

test_that("frmtmb_register_aterm validates its own arguments", {
  expect_error(frmtmb_register_aterm("trials"),
               "cannot be re-registered")
  expect_error(frmtmb_register_aterm("has space"), "one syntactic name")
  expect_error(frmtmb_register_aterm("ok", arity = 0),
               "whole number of arguments")
  expect_error(frmtmb_register_aterm("ok", coerce = "as.numeric"),
               "must be a function of one")
})

test_that("a registered coercion must return numbers", {
  on.exit(frmtmb_aterm_registry$reg[["tag"]] <- NULL, add = TRUE)
  frmtmb_register_aterm("tag", coerce = function(x) as.character(x))
  set.seed(33)
  dd <- data.frame(x = stats::rnorm(30), flag = rep(0:1, 15))
  dd$y <- stats::rnorm(30, dd$x)
  expect_error(frm(bf(y | tag(flag) ~ x) + gaussian(), data = dd),
               "must be numeric")
})
