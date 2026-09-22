# brms's return shapes on DRAWS (item 2.6f): the three summarizing
# methods, the two deprecated accessors brms keeps live, and
# point_estimate.
#
# Run against frmtmb.sample 0.8.0 first and recorded failing in
# dev/shapes-log/seen-failing-sample.txt: fitted() and residuals()
# returned NULL through stats::fitted.default(), nsamples() and
# posterior_samples() were refused, and point_estimate was ignored.

shapes_draws <- local({
  cache <- new.env(parent = emptyenv())
  function() {
    if (!is.null(cache$ds)) return(cache$ds)
    set.seed(20260917)
    n <- 60
    dd <- data.frame(x = rnorm(n))
    dd$y <- rnorm(n, 1 + 0.5 * dd$x, 1)
    fit <- frm(bf(y ~ x) + gaussian(), data = dd)
    cache$ds <- frm_sample(fit, chains = 1, iter = 150, warmup = 100,
                           seed = 20260917, refresh = 0)
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
  # frm_sample() stores no warmup, so there is nothing to count there
  expect_error(suppressWarnings(nsamples(ds, incl_warmup = TRUE)),
               "warmup")
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
  ds <- frm_sample(fit, chains = 1, iter = 300, warmup = 150,
                   seed = 20260921, refresh = 0)
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
  ds <- frm_sample(fit, chains = 1, iter = 200, warmup = 100,
                   seed = 20260922, refresh = 0)
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
