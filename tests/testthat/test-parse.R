test_that("bf() builds a formula object and + attaches a family", {
  f <- bf(y ~ x + (x | g))
  expect_s3_class(f, "frmtmb_formula")
  expect_null(f$family)
  f2 <- f + gaussian()
  expect_s3_class(f2$family, "frmtmb_family")
  expect_identical(f2$family$family, "gaussian")
})

test_that("bf() rejects unsupported features with clear errors", {
  # nl = TRUE without a parameter formula waits for frm(), as brms
  # waits for brm()
  expect_error(frm(bf(y ~ a * exp(b * x), nl = TRUE), data = NULL,
                   dry_run = "spec"),
               "nonlinear-parameter formula")
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
  # one term, as brms reads (x | gr(g, cor = FALSE))
  spec <- frm(bf(y ~ x + (x || g)) + poisson(),
                 data = NULL, dry_run = "spec")
  re <- spec$responses[[1]]$dpars$mu$re
  expect_length(re, 1)
  expect_identical(re[[1]]$covstruct, "diag")

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
  expect_error(frm(bf(y | side(d) ~ x) + gaussian(),
                   data = NULL, dry_run = "spec"),
               "side")
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

# --- the addition-term registry (frmtmb.eam finding 1) ----------------
# `side()` and not `dec()`: frmtmb.eam registers dec() when it loads,
# and covr runs every file in one session, so an earlier file that
# loads frmtmb.eam made dec() parse here and the on.exit() that was
# here removed frmtmb.eam's entry for the files after this one. Each
# test now saves the registries and restores them (helper-registry.R):
# a deletion by name would also leave the side() compat feature behind,
# since frmtmb_register_aterm() writes one.

test_that("the refusal names vint()/vreal() and the registration seam", {
  # an extension author reading this message has to learn that the set
  # is not closed and what to reach for meanwhile
  err <- tryCatch(frm(bf(y | side(bound) ~ x) + gaussian(), data = NULL,
                      dry_run = "spec"),
                  error = conditionMessage)
  expect_match(err, "vint() for integers", fixed = TRUE)
  expect_match(err, "frmtmb_register_aterm()", fixed = TRUE)
})

test_that("a registered aterm parses, coerces and reaches the density", {
  local_frmtmb_registries()
  frmtmb_register_aterm("side", arity = 1,
                        coerce = function(x) as.integer(factor(x)) - 1L)

  spec <- frm(bf(y | side(bound) ~ x) + gaussian(), data = NULL,
              dry_run = "spec")
  expect_identical(deparse1(spec$responses[[1]]$aterms$side), "bound")

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
    "side_reader", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu + 0.8 * aterms$side, dpars$sigma,
                  log = TRUE)
    },
    required_aterms = "side",
    post = list(mean_fn = function(dpars, aterms) {
      dpars$mu + 0.8 * aterms$side
    }))
  fit <- frm(bf(y | side(bound) ~ x) + fam, data = dd)
  # the factor arrived as the 0/1 the coercion produced, not as level
  # codes and not as a factor
  seen <- fit$frame$aterm_values$y$side
  expect_identical(sort(unique(seen)), c(0, 1))
  expect_identical(seen[1:2], c(0, 1))
  # and the same coercion runs on newdata
  nd <- data.frame(x = c(0, 0),
                   bound = factor(c("upper", "lower"),
                                  levels = levels(dd$bound)))
  pr <- frm_linpred(fit, newdata = nd, type = "response")
  expect_gt(pr[1], pr[2])
})

test_that("a registered aterm is required on newdata that omits it", {
  local_frmtmb_registries()
  frmtmb_register_aterm("side")
  set.seed(32)
  dd <- data.frame(x = stats::rnorm(60), bound = rep(0:1, 30))
  dd$y <- stats::rnorm(60, dd$x + dd$bound)
  fam <- custom_family(
    "side_reader2", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu + aterms$side, dpars$sigma, log = TRUE)
    },
    post = list(mean_fn = function(dpars, aterms) dpars$mu + aterms$side))
  fit <- frm(bf(y | side(bound) ~ x) + fam, data = dd)
  expect_error(frm_linpred(fit, newdata = data.frame(x = 0),
                       type = "response"),
               "could not be evaluated on newdata")
})

test_that("a registered aterm of arity above one numbers its arguments", {
  local_frmtmb_registries()
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

test_that("a registration whose key another term delivers under is refused", {
  local_frmtmb_registries()
  frmtmb_register_aterm("tterm", arity = 3)
  reg0 <- frmtmb_aterm_registry$reg
  feat0 <- frm_compat_features()
  # tterm2 is the key of tterm()'s second argument: through 0.81.0 this
  # was accepted, and `y | tterm2(z) + tterm(u, v, w)` handed v to the
  # family as tterm2
  expect_error(frmtmb_register_aterm("tterm2"),
               paste0("under the key `tterm2`, which the registered ",
                      "`tterm()` (3 arguments) already uses"),
               fixed = TRUE)
  # the refusal left both registries as it found them
  expect_identical(frmtmb_aterm_registry$reg, reg0)
  expect_identical(frm_compat_features(), feat0)
  # the keys core delivers its own terms' values under
  expect_error(frmtmb_register_aterm("vint1"),
               "which frmtmb's own `vint()` already uses", fixed = TRUE)
  expect_error(frmtmb_register_aterm("cens_y", arity = 2),
               "key `cens_y2`, which frmtmb's own `cens()`", fixed = TRUE)
  expect_error(frmtmb_register_aterm("trunc_lb"),
               "which frmtmb's own `trunc()` already uses", fixed = TRUE)
  # the absent case: names that share a stem but no key, and the same
  # name again, which replaces it
  expect_no_error(frmtmb_register_aterm("tterm1x"))
  expect_no_error(frmtmb_register_aterm("ttermb", arity = 2))
  expect_no_error(frmtmb_register_aterm("tterm", arity = 3))
  spec <- frm(bf(y | ttermb(a, b) + tterm(u, v, w) ~ x) + gaussian(),
              data = NULL, dry_run = "spec")
  expect_named(spec$responses[[1]]$aterms,
               c("ttermb1", "ttermb2", "tterm1", "tterm2", "tterm3"))
})

test_that("a registration whose keys reach a registered name is refused", {
  local_frmtmb_registries()
  # the other order: the earlier term owns the key `tterm2`. Through
  # 0.81.0 `y | tterm2(z) + tterm(u, v, w)` then stopped with
  # "Duplicated addition term `tterm2()`", which names neither cause
  frmtmb_register_aterm("tterm2")
  expect_error(frmtmb_register_aterm("tterm", arity = 3),
               paste0("`tterm()` with 3 arguments would deliver a value ",
                      "under the key `tterm2`"),
               fixed = TRUE)
  # twelve arguments reach ta11 and ta12, a two-argument ta1()'s keys
  frmtmb_register_aterm("ta1", arity = 2)
  expect_error(frmtmb_register_aterm("ta", arity = 12),
               "key `ta11`, which the registered `ta1()` (2 arguments)",
               fixed = TRUE)
  # an option term is found by its name, so its name is a key
  frmtmb_register_aterm("tq2", option = function(family, call, env) {
    family
  })
  expect_error(frmtmb_register_aterm("tq", arity = 2),
               "key `tq2`, which the registered `tq2()` already uses",
               fixed = TRUE)
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
  local_frmtmb_registries()
  frmtmb_register_aterm("tag", coerce = function(x) as.character(x))
  set.seed(33)
  dd <- data.frame(x = stats::rnorm(30), flag = rep(0:1, 15))
  dd$y <- stats::rnorm(30, dd$x)
  expect_error(frm(bf(y | tag(flag) ~ x) + gaussian(), data = dd),
               "must be numeric")
})
