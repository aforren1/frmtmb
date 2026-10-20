# brms's sd_ and cor_ variables on frm_sample()'s draws. The sampler
# works on the internal covariance parameters (theta_1, ...), and
# through frmtmb.sample 0.19.0 those were the only covariance columns,
# so brms-style code such as as_draws_df(draws)$sd_g__Intercept found
# nothing. The derived columns take their names from core's
# hyp_env_vals(), which is also where variables() of a fit takes them.

skip_on_cran()

sdn_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    skip_if_not_installed("posterior")
    if (is.null(cache)) {
      set.seed(11)
      n <- 160
      dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:16, 10)),
                       h = factor(rep(1:8, each = 20)))
      u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
      dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + u[dd$g, 1] +
                             u[dd$g, 2] * dd$x +
                             stats::rnorm(8, 0, 0.5)[dd$h], 1)
      fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), family = gaussian(),
                 data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 400, refresh = 0, seed = 3)))
      cache <<- list(fit = fit, ds = ds)
    }
    cache
  }
})

test_that("the draws carry brms's sd_ and cor_ variables", {
  cs <- sdn_case()
  v <- posterior::variables(cs$ds)
  want <- grep("^(sd|cor)_", variables(cs$fit), value = TRUE)
  # brms 2.23.0's names for this model, in its order: every sd_, then
  # every cor_ (lane prilabel, dev/prilabel-findings.md)
  expect_identical(want, c("sd_g__Intercept", "sd_g__x",
                           "sd_h__Intercept", "cor_g__Intercept__x"))
  expect_true(all(want %in% v))
  adf <- posterior::as_draws_df(cs$ds)
  expect_true(all(want %in% names(adf)))
  # the values are the draws' own: exp(theta) for a standard deviation,
  # and VarCorr()'s for both, draw by draw
  m <- as.matrix(cs$ds)
  expect_identical(adf$sd_h__Intercept, unname(exp(m[, "theta_4"])))
  vc <- VarCorr(cs$ds, summary = FALSE)
  expect_equal(adf$sd_g__x, unname(vc$g$sd[, "x"]), tolerance = 1e-14)
  expect_equal(adf$cor_g__Intercept__x, unname(vc$g$cor[, "Intercept", "x"]),
               tolerance = 1e-14)
})

test_that("the derived columns stay out of what the model reads", {
  cs <- sdn_case()
  im <- frmtmb.sample:::draws_internal_matrix(cs$ds)
  expect_false(any(grepl("^(sd|cor)_", colnames(im))))
  # every method that hands a draw back to the model still runs, and
  # summary() reports the derived rows beside the internal ones
  ep <- posterior_epred(cs$ds, ndraws = 5)
  expect_identical(dim(ep), c(5L, 160L))
  s <- summary(cs$ds)
  expect_true(all(c("sd_g__Intercept", "cor_g__Intercept__x") %in%
                    rownames(s)))
})

test_that("the fast path reproduces hyp_env_vals() or steps aside", {
  cs <- sdn_case()
  nc <- frmtmb.sample:::draws_natural_cols(cs$ds$fit)
  expect_false(is.null(nc$derived$fast))
  expect_setequal(nc$derived$fast$names, nc$derived$names)
  # the slow path, every draw through hyp_env_vals(), gives the same
  # columns: the guard's absent case
  slow <- nc$derived
  slow$fast <- NULL
  im <- frmtmb.sample:::draws_internal_matrix(cs$ds, rows = 1:20)
  a <- frmtmb.sample:::draws_cov_values(im, cs$ds$fit, nc$derived)
  b <- frmtmb.sample:::draws_cov_values(im, cs$ds$fit, slow)
  expect_equal(a, b, tolerance = 1e-14)
})
