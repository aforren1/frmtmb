# brms's return shapes on DRAWS (item 2.6f): the three summarizing
# methods, the two deprecated accessors brms keeps live, and
# point_estimate.
#
# Run against frmtmb.sample 0.8.0 first and recorded failing in
# dev/shapes-log/seen-failing-sample.txt: fitted() and residuals()
# returned NULL through stats::fitted.default(), nsamples() and
# posterior_samples() were refused, and point_estimate was ignored.

# The chains here are short on purpose: these tests are about the SHAPE
# of what the draws methods return, and rstan says the effective sample
# size is low, which is true and beside the point.
short_chain <- c("Bulk Effective Samples Size",
                 "Tail Effective Samples Size")

shapes_draws <- local({
  cache <- new.env(parent = emptyenv())
  function() {
    if (!is.null(cache$ds)) return(cache$ds)
    set.seed(20260917)
    n <- 60
    dd <- data.frame(x = rnorm(n))
    dd$y <- rnorm(n, 1 + 0.5 * dd$x, 1)
    fit <- frm(bf(y ~ x) + gaussian(), data = dd)
    cache$ds <- allow_warnings(
      frm_sample(fit, chains = 1, iter = 150, warmup = 100,
                 seed = 20260917, refresh = 0),
      short_chain)
    cache$ds
  }
})

test_that("fitted() on draws summarizes posterior_epred()", {
  skip_on_cran()
  skip_sampler()
  ds <- shapes_draws()
  fi <- fitted(ds)
  expect_equal(dim(fi), c(nobs(ds$fit), 4L))
  expect_equal(colnames(fi), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  ep <- posterior_epred(ds)
  expect_equal(unname(fi[, "Estimate"]), unname(colMeans(ep)))
  expect_equal(dim(fitted(ds, summary = FALSE)), dim(ep))
})

test_that("predict() and residuals() on draws are brms's summaries", {
  skip_on_cran()
  skip_sampler()
  ds <- shapes_draws()
  set.seed(1)
  pr <- predict(ds)
  expect_equal(dim(pr), c(nobs(ds$fit), 4L))
  expect_equal(colnames(pr), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  rs <- residuals(ds)
  expect_equal(dim(rs), c(nobs(ds$fit), 4L))
  # the predictive interval is far wider than the expected-value one
  expect_gt(mean(pr[, "Q97.5"] - pr[, "Q2.5"]),
            3 * mean(fitted(ds)[, "Q97.5"] - fitted(ds)[, "Q2.5"]))
})

test_that("point_estimate collapses the draws, as in brms", {
  skip_on_cran()
  skip_sampler()
  ds <- shapes_draws()
  ep <- posterior_epred(ds, point_estimate = "median",
                        ndraws_point_estimate = 2)
  expect_equal(nrow(ep), 2L)
  # the two rows are the SAME parameter vector, so they agree exactly
  expect_identical(ep[1, ], ep[2, ])
  expect_equal(nrow(posterior_epred(ds, point_estimate = "mean")), 1L)
  expect_equal(nrow(posterior_predict(ds, point_estimate = "mean")), 1L)
  expect_equal(nrow(posterior_linpred(ds, point_estimate = "mean")), 1L)
  expect_equal(nrow(log_lik(ds, point_estimate = "mean")), 1L)
  expect_error(posterior_epred(ds, point_estimate = "mode"), "mode")
})

test_that("nsamples() is brms's, with brms's deprecation warning", {
  skip_on_cran()
  skip_sampler()
  ds <- shapes_draws()
  expect_warning(n <- nsamples(ds), "deprecated")
  expect_equal(n, ndraws(ds))
  expect_equal(suppressWarnings(nsamples(ds, subset = 10:1)), 10L)
  # brmsfit-methods:595: the saved iterations, warmup included, from the
  # stanfit's record, as brms counts them. 0.66.0's frmtmb.sample refused
  # this ("has nothing to count"). shapes_draws() runs one chain of
  # iter = 150 with warmup = 100, and rstan saves the warmup
  expect_equal(suppressWarnings(nsamples(ds, incl_warmup = TRUE)), 150L)
  expect_identical(ndraws(ds), 50L)
  # brms's subset check reads the count asked for
  expect_error(suppressWarnings(nsamples(ds, subset = 1:60)),
               "Argument 'subset' is invalid.", fixed = TRUE)
  expect_equal(suppressWarnings(nsamples(ds, subset = 1:60,
                                         incl_warmup = TRUE)), 60L)
  # draws with no stanfit behind them have no record to read
  bare <- ds
  bare$stanfit <- NULL
  expect_error(suppressWarnings(nsamples(bare, incl_warmup = TRUE)),
               "carry no stanfit")
  expect_equal(suppressWarnings(nsamples(bare)), ndraws(ds))
})

test_that("posterior_samples() is brms's, with brms's warning", {
  skip_on_cran()
  skip_sampler()
  ds <- shapes_draws()
  expect_warning(d <- posterior_samples(ds), "deprecated")
  expect_true(is.data.frame(d))
  expect_equal(dim(d), c(ndraws(ds), length(variables(ds))))
  expect_equal(names(d), variables(ds))
  # `pars` is a regular expression unless fixed = TRUE, as in brms
  b <- suppressWarnings(posterior_samples(ds, pars = "^b_"))
  expect_equal(names(b), grep("^b_", variables(ds), value = TRUE))
  expect_true(is.matrix(suppressWarnings(
    posterior_samples(ds, as.matrix = TRUE))))
})

test_that("posterior_samples(pars = ) orders the coefficients as brms", {
  # brmsfit-methods:635: brms lists every intercept first. 0.66.0 gave
  # variables() order, each predictor's coefficients together:
  # b_Intercept, b_x, b_sigma_Intercept, b_sigma_x
  set.seed(20260930)
  n <- 80
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x, exp(0.2 + 0.3 * dd$x))
  fit <- frm(bf(y ~ x, sigma ~ x), family = gaussian(), data = dd)
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  M <- matrix(rep(est, each = 6) + stats::rnorm(6 * length(est), 0, 0.05),
              6, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  ds <- structure(list(stanfit = NULL,
                       draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                                     lp__ = 0),
                       fit = fit), class = "frmtmb_draws")
  ps <- function(...) {
    allow_warnings(posterior_samples(ds, ...), "is deprecated")
  }
  brms_b <- c("b_Intercept", "b_sigma_Intercept", "b_x", "b_sigma_x")
  expect_identical(names(ps(pars = "^b_")), brms_b)
  expect_identical(paste0("b_", rownames(fixef(fit))), brms_b)
  # the values move with their names
  expect_identical(ps(pars = "^b_")$b_sigma_Intercept,
                   unname(ds$draws[, "b_sigma_Intercept"]))
  # several patterns: each keeps brms's order within its matches
  expect_identical(names(ps(pars = c("sigma", "^b_x$"))),
                   c("b_sigma_Intercept", "b_sigma_x", "b_x"))
  # fixed = TRUE keeps the order given, as brms's intersect() does
  expect_identical(names(ps(pars = c("b_x", "b_Intercept"), fixed = TRUE)),
                   c("b_x", "b_Intercept"))
  # the guard absent: no pars is every variable, in variables() order
  expect_identical(names(ps()), variables(ds))
  # and a pattern that matches nothing is brms's NULL
  expect_null(ps(pars = "^nothing_"))
})

test_that("parnames() is variables(), with brms's warning", {
  skip_on_cran()
  skip_sampler()
  ds <- shapes_draws()
  expect_warning(p <- parnames(ds), "deprecated")
  expect_equal(p, variables(ds))
})

# Punch round 2.

test_that("fixef() on ordinal draws has the thresholds, in brms's order", {
  skip_on_cran()
  skip_sampler()
  set.seed(20260921)
  n <- 120
  dd <- data.frame(x = rnorm(n))
  dd$ord <- factor(cut(0.8 * dd$x + stats::rlogis(n),
                       c(-Inf, -0.5, 0.8, Inf), labels = 1:3),
                   ordered = TRUE)
  fit <- frm(bf(ord ~ x) + cumulative(), data = dd)
  ds <- allow_warnings(frm_sample(fit, chains = 1, iter = 300, warmup = 150,
                                  seed = 20260921, refresh = 0),
                       short_chain)
  fe <- fixef(ds)
  # the fit's rows and brms's: it used to report `x` alone
  expect_equal(rownames(fe), rownames(fixef(fit)))
  expect_equal(rownames(fe), c("Intercept[1]", "Intercept[2]", "x"))
  # the threshold is the model's own, not the internal log increment
  expect_gt(fe["Intercept[2]", "Estimate"], fe["Intercept[1]", "Estimate"])
  d <- fixef(ds, summary = FALSE)
  expect_true(all(d[, "Intercept[2]"] > d[, "Intercept[1]"]))
})

test_that("allow_new_levels on draws refuses only an unseen level", {
  skip_on_cran()
  skip_sampler()
  set.seed(20260922)
  n <- 60
  dd <- data.frame(x = rnorm(n), g = factor(rep(1:6, each = 10)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(6, 0, 0.7)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- allow_warnings(frm_sample(fit, chains = 1, iter = 200, warmup = 100,
                                  seed = 20260922, refresh = 0),
                       short_chain)
  known <- data.frame(x = c(0, 1), g = factor(c("1", "2")))
  unseen <- data.frame(x = 0, g = factor("new"))
  for (fn in c("posterior_predict", "posterior_epred", "posterior_linpred",
               "predict", "fitted")) {
    run <- function(...) {
      set.seed(1)
      do.call(fn, list(ds, ...))
    }
    # brms's default, FALSE, answers as it did before
    expect_equal(run(allow_new_levels = FALSE), run(), info = fn)
    # TRUE with no newdata, or with levels the fit saw, changes nothing
    expect_equal(run(allow_new_levels = TRUE), run(), info = fn)
    expect_equal(run(newdata = known, allow_new_levels = TRUE),
                 run(newdata = known), info = fn)
    # TRUE with an unseen level is refused, naming the function CALLED
    err <- tryCatch(run(newdata = unseen, allow_new_levels = TRUE),
                    error = function(e) conditionMessage(e))
    expect_true(is.character(err), info = fn)
    expect_match(err, paste0("^", fn, "[(][)] on draws"), info = fn)
    expect_match(err, "neither is implemented for frmtmb_draws", info = fn)
  }
})

test_that("posterior_samples(pars = ) orders a nonlinear model as brms", {
  # review m5: fixef()'s order put every intercept first, a nonlinear
  # parameter's too. brms 2.23.0 keeps a nonlinear parameter's intercept
  # with its coefficients and puts a distributional parameter's first
  # (dev/ceplot-log/p1-psorder-brms.txt)
  set.seed(11)
  n <- 200
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
  d$yp <- 2 * exp(0.3 * d$x) + 0.2 * d$z + rnorm(n, 0, exp(0.1 * d$w))
  hand <- function(fit) {
    tpl <- fit$frame$par_template
    est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
    M <- matrix(rep(est, each = 4) + stats::rnorm(4 * length(est), 0, 0.01),
                4, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
    structure(list(stanfit = NULL,
                   draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                                 lp__ = 0),
                   fit = fit), class = "frmtmb_draws")
  }
  ps <- function(ds, ...) {
    allow_warnings(names(posterior_samples(ds, ...)), "is deprecated")
  }
  f1 <- frm(bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
            family = gaussian(), data = d)
  expect_identical(ps(hand(f1), pars = "^b_"),
                   c("b_a_Intercept", "b_a_z", "b_b_Intercept"))
  f2 <- frm(bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1 + w, sigma ~ w,
               nl = TRUE), family = gaussian(), data = d)
  expect_identical(ps(hand(f2), pars = "^b_"),
                   c("b_sigma_Intercept", "b_a_Intercept", "b_a_z",
                     "b_b_Intercept", "b_b_w", "b_sigma_w"))
  f3 <- frm(bf(yp ~ x + z, sigma ~ w), family = gaussian(), data = d)
  expect_identical(ps(hand(f3), pars = "^b_"),
                   c("b_Intercept", "b_sigma_Intercept", "b_x", "b_z",
                     "b_sigma_w"))
})
