# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived, and this file is a main contributor
# to CRAN-condition check time.
skip_on_cran()

# v0.7: OSA residuals, interval censoring, discrete truncation,
# sampling bridge, natural-scale varcorr CIs, smooth edf.

test_that("OSA residuals are standard normal for correct models", {
  data(sleepstudy, package = "lme4")
  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
             data = sleepstudy)
  r <- residuals(fit, type = "osa")
  expect_length(r, 180)
  expect_gt(stats::ks.test(r, "pnorm")$p.value, 0.01)

  set.seed(201)
  dp <- data.frame(x = rnorm(300), g = factor(rep(1:15, 20)))
  dp$y <- rpois(300, exp(0.4 + 0.3 * dp$x + rnorm(15, 0, 0.4)[dp$g]))
  fitp <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dp)
  rp <- residuals(fitp, type = "osa")
  expect_gt(stats::ks.test(rp, "pnorm")$p.value, 0.01)
})

test_that("interval censoring matches a hand-rolled reference (brms#1070)", {
  set.seed(202)
  n <- 400
  x <- rnorm(n)
  ystar <- 1 + 0.6 * x + rnorm(n)
  # interval-censor a third of the observations into unit bins
  cen <- ifelse(seq_len(n) %% 3 == 0, 2, 0)
  ylo <- ifelse(cen == 2, floor(ystar), ystar)
  yhi <- ifelse(cen == 2, floor(ystar) + 1, NA)   # NA on non-interval rows
  dd <- data.frame(y = ylo, y2 = yhi, x = x, cen = cen)

  fit <- frm(bf(y | cens(cen, y2) ~ x) + gaussian(), data = dd)

  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * dd$x
    s <- exp(p$ls)
    ll <- RTMB::dnorm(dd$y, mu, s, log = TRUE)
    i2 <- which(dd$cen == 2)
    ll[i2] <- log(RTMB::pnorm((dd$y[i2] + 1 - mu[i2]) / s) -
                    RTMB::pnorm((dd$y[i2] - mu[i2]) / s))
    -sum(ll)
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)

  # validations
  dd_bad <- dd; dd_bad$y2[dd_bad$cen == 2][1] <- NA
  expect_error(frm(bf(y | cens(cen, y2) ~ x) + gaussian(), data = dd_bad),
               "must not be NA")
  expect_error(frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd),
               "needs upper bounds")
})

test_that("zero-truncated poisson matches glmmTMB truncated_poisson", {
  skip_if_not_installed("glmmTMB")
  set.seed(203)
  n <- 2000
  x <- rnorm(n)
  y <- rpois(n, exp(0.6 + 0.4 * x))
  keep <- y > 0
  dd <- data.frame(y = y[keep], x = x[keep])
  fit <- suppressWarnings(
    frm(bf(y | trunc(lb = 1) ~ x) + poisson(), data = dd)
  )
  ref <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::truncated_poisson,
                          data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
  expect_error(frm(bf(y | trunc(lb = 0) ~ x) + poisson(), data = dd),
               "lb >= 1")
})

test_that("frm_sample returns named draws and check_laplace agrees on a clean model", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(204)
  dd <- data.frame(x = rnorm(150), g = factor(rep(1:15, 10)))
  dd$y <- rnorm(150, 1 + 0.5 * dd$x + rnorm(15, 0, 0.8)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 600,
                                    refresh = 0, seed = 1))
  m <- as.matrix(ds)
  # draws names are parenthesis-free (the brms convention; v0.36)
  expect_true(all(c("Intercept", "x", "sigma_Intercept",
                    "theta_1") %in% colnames(m)))
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and this asserts wiring, not mixing
  if (sampler_gates_on()) {
    expect_lt(abs(mean(m[, "x"]) - fixef(fit)$mu[["x"]]),
              5 * stats::sd(m[, "x"]) + 1e-8)
  }

  cl <- suppressWarnings(suppressMessages(
    check_laplace(fit, chains = 1, iter = 600, refresh = 0, seed = 1)))
  expect_s3_class(cl, "data.frame")
  expect_true("ess_bulk" %in% names(cl))
  # Wald and posterior agree on a well-behaved gaussian LMM, but only a
  # HEALTHY chain can testify: on a platform whose chain wandered
  # (measured, not assumed), the agreement claim is untestable
  row_x <- cl[cl$parameter == "x", ]
  # bulk ESS is necessary, not sufficient: a chain can mix on x while
  # its flat-prior theta excursion fattens the marginal anyway, so the
  # agreement claim is additionally gated per platform
  if (sampler_gates_on() &&
      is.finite(row_x$ess_bulk) && row_x$ess_bulk >= 100) {
    expect_lt(abs(row_x$z_shift), 0.75)
    expect_lt(abs(row_x$sd_ratio - 1), 0.5)
  } else {
    skip("chain too unhealthy on this platform to judge the agreement")
  }
})

test_that("confint_varcorr matches glmmTMB's natural-scale intervals", {
  skip_if_not_installed("glmmTMB")
  data(sleepstudy, package = "lme4")
  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
             data = sleepstudy)
  cv <- confint_varcorr(fit)
  expect_identical(cv$type, c("sd", "sd", "cor"))
  expect_true(all(cv$lwr < cv$estimate & cv$estimate < cv$upr))

  ref <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject),
                          sleepstudy, REML = FALSE)
  ci_g <- confint(ref)   # wald, natural scale for sds
  sd_rows <- grep("Std\\.Dev", rownames(ci_g))
  expect_vector_equal(cv$estimate[1:2], unname(ci_g[sd_rows, "Estimate"]),
                      tol = 1e-3)
  expect_vector_equal(cv$lwr[1:2], unname(ci_g[sd_rows, "2.5 %"]),
                      tol = 0.05)
  expect_vector_equal(cv$upr[1:2], unname(ci_g[sd_rows, "97.5 %"]),
                      tol = 0.05)
})

test_that("smooth edf approximates mgcv's", {
  set.seed(205)
  n <- 300
  dd <- data.frame(x = runif(n))
  dd$y <- sin(3 * dd$x) + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ s(x)) + gaussian(), data = dd)
  s <- summary(fit)
  expect_false(is.null(s$smooth_edf))
  ref <- mgcv::gam(y ~ s(x), data = dd, method = "ML")
  edf_gam <- sum(ref$edf) - 1   # gam's s(x) edf (intercept removed)
  # ours counts the penalized part; add the 1 null-space column to align
  expect_lt(abs((s$smooth_edf[[1]] + 1) - edf_gam), 1)
})

test_that("marginaleffects works through the extension generics", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("glmmTMB")
  set.seed(206)
  dd <- data.frame(x = rnorm(300),
                   f = factor(sample(c("a", "b"), 300, replace = TRUE)),
                   g = factor(rep(1:15, 20)))
  dd$y <- rpois(300, exp(0.3 + 0.4 * dd$x + 0.5 * (dd$f == "b") +
                           rnorm(15, 0, 0.4)[dd$g]))
  fit <- frm(bf(y ~ x + f + (1 | g)) + poisson(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x + f + (1 | g), family = poisson,
                          data = dd)

  # variables= keeps marginaleffects off the grouping factor (it has no
  # way to know g is a random-effect group for an external class)
  sl_f <- marginaleffects::avg_slopes(fit, newdata = dd,
                                      variables = c("x", "f"))
  sl_g <- suppressWarnings(
    marginaleffects::avg_slopes(ref, newdata = dd,
                                variables = c("x", "f"))
  )
  of <- order(sl_f$term); og <- order(sl_g$term)
  expect_vector_equal(sl_f$estimate[of], sl_g$estimate[og], tol = 1e-4)
  expect_vector_equal(sl_f$std.error[of], sl_g$std.error[og], tol = 1e-3)

  cm_f <- marginaleffects::avg_comparisons(fit, newdata = dd,
                                           variables = "f")
  expect_true(is.finite(cm_f$estimate))
})