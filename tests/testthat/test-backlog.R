# Second wave of mined edge cases (dev/test-backlog.md).

test_that("non-default contrasts survive prediction (glmmTMB#439)", {
  set.seed(161)
  dd <- data.frame(f = factor(rep(c("a", "b", "c"), 40)),
                   g = factor(rep(1:6, 20)))
  dd$y <- rnorm(120, c(a = 0, b = 1, c = -1)[dd$f] +
                  rnorm(6, 0, 0.5)[dd$g], 1)
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  fit <- frm(bf(y ~ f + (1 | g)) + gaussian(), data = dd)
  p0 <- predict(fit)
  options(old)   # flip back BEFORE predicting on newdata
  expect_equal(predict(fit, newdata = dd), p0, tolerance = 1e-8)
  # explicit contrasts= attribute on the factor
  dd2 <- dd
  contrasts(dd2$f) <- stats::contr.sum(3)
  fit2 <- frm(bf(y ~ f + (1 | g)) + gaussian(), data = dd2)
  # xlev refactoring drops the factor's contrasts attr with a warning;
  # the stored contrasts matrix still applies, so values must match
  expect_equal(suppressWarnings(predict(fit2, newdata = dd2)),
               predict(fit2), tolerance = 1e-8)
})

test_that("re.form = NA needs no grouping columns in newdata (glmmTMB#923)", {
  set.seed(162)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.6)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  nd <- data.frame(x = c(-1, 0, 1))    # no g column at all
  p <- predict(fit, newdata = nd, re.form = NA, se.fit = TRUE)
  expect_length(p$fit, 3)
  expect_true(all(is.finite(p$se.fit)))
})

test_that("fits survive the calling environment disappearing (lme4 formulaEval)", {
  set.seed(163)
  make_fit <- function() {
    local_dat <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
    local_dat$y <- rnorm(80, 1 + 0.5 * local_dat$x, 1)
    fml <- stats::as.formula(paste("y ~", "x + (1 | g)"))
    frm(bf(fml) + gaussian(), data = local_dat)
  }
  fit <- make_fit()
  # local_dat is gone; the stored frame must carry prediction and emmeans
  expect_length(predict(fit, newdata = model.frame(fit)), 80)
  expect_length(fitted(fit), 80)
  skip_if_not_installed("emmeans")
  em <- as.data.frame(emmeans::emmeans(fit, ~1))
  expect_true(is.finite(em$emmean[1]))
})

test_that("slash and interaction grouping syntax expand correctly", {
  set.seed(164)
  dd <- expand.grid(a = factor(1:6), b = factor(1:4), rep = 1:8)
  dd$y <- rnorm(nrow(dd),
                rnorm(6, 0, 0.6)[dd$a] +
                  rnorm(24, 0, 0.4)[as.integer(interaction(dd$a, dd$b))],
                1)
  f1 <- frm(bf(y ~ 1 + (1 | a / b)) + gaussian(), data = dd)
  f2 <- frm(bf(y ~ 1 + (1 | a) + (1 | a:b)) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-6)
  expect_length(f1$frame$re_blocks, 2)
  # the interaction block carries one level per observed combination
  bk <- Filter(function(b) grepl(":", b$group_name),
               f1$frame$re_blocks)[[1]]
  expect_length(bk$levels, nlevels(droplevels(dd$a:dd$b)))
})

test_that("cens() composes with dpar formulas (brms#1716)", {
  set.seed(165)
  n <- 300
  x <- rnorm(n); z <- rnorm(n)
  ystar <- 1 + 0.5 * x + rnorm(n, 0, exp(0.2 + 0.3 * z))
  cp <- 2
  dd <- data.frame(y = pmin(ystar, cp), x = x, z = z,
                   cen = as.numeric(ystar > cp))
  fit <- frm(bf(y | cens(cen) ~ x, sigma ~ z) + gaussian(), data = dd)

  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * dd$x
    s <- exp(p$bs[1] + p$bs[2] * dd$z)
    ll <- RTMB::dnorm(dd$y, mu, s, log = TRUE)
    ir <- which(dd$cen == 1)
    ll[ir] <- log(1 - RTMB::pnorm((dd$y[ir] - mu[ir]) / s[ir]))
    -sum(ll)
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), bs = c(0, 0)),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("weights are frequency-like: aggregated poisson equivalence", {
  set.seed(166)
  long <- data.frame(x = rep(c(0, 1), each = 200))
  long$y <- rpois(400, exp(0.5 + 0.4 * long$x))
  agg <- aggregate(cnt ~ y + x, transform(long, cnt = 1), sum)
  m1 <- frm(bf(y ~ x) + poisson(), data = long)
  m2 <- frm(bf(y | weights(cnt) ~ x) + poisson(), data = agg)
  expect_vector_equal(fixef(m1)$mu, fixef(m2)$mu, tol = 1e-5)
  expect_lt(abs(as.numeric(logLik(m1)) - as.numeric(logLik(m2))), 1e-5)
})

test_that("Inf responses error; dpar names with underscores rejected", {
  dd <- data.frame(y = c(1, Inf, 3), x = 1:3)
  expect_error(frm(bf(y ~ x) + gaussian(), data = dd), "Non-finite")
  expect_error(bf(y ~ x, my_par ~ z), "dots or underscores")
})

test_that("predict warns on unknown arguments", {
  dd <- data.frame(y = rnorm(30), x = rnorm(30))
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  expect_warning(predict(fit, bogus = 1), "ignoring unknown arguments")
  # the common typo se= partial-matches se.fit and just works
  p <- predict(fit, se = TRUE)
  expect_named(p, c("fit", "se.fit"))
})

test_that("na.exclude pads fitted/residuals/predict to full length", {
  set.seed(167)
  dd <- data.frame(x = rnorm(50), g = factor(rep(1:5, 10)))
  dd$y <- rnorm(50, 1 + 0.5 * dd$x, 1)
  dd$y[c(3, 7)] <- NA
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
             na.action = stats::na.exclude)
  expect_identical(stats::nobs(fit), 48L)
  expect_length(fitted(fit), 50)
  expect_true(all(is.na(fitted(fit)[c(3, 7)])))
  expect_length(residuals(fit), 50)
  expect_length(predict(fit), 50)
  # na.omit(padded) equals the na.omit fit's values
  fit0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_equal(unname(as.vector(stats::na.omit(fitted(fit)))),
               unname(fitted(fit0)), tolerance = 1e-8)
})

test_that("lazy sdreport: estimates identical, SEs on demand", {
  set.seed(168)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.6)[dd$g], 1)
  f_lazy <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  f_eager <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, se = TRUE)
  expect_null(f_lazy$cache$sdr)
  expect_false(is.null(f_eager$cache$sdr))
  expect_vector_equal(fixef(f_lazy)$mu, fixef(f_eager)$mu, tol = 1e-10)
  expect_vector_equal(unlist(f_lazy$estimates$b),
                      unlist(f_eager$estimates$b), tol = 1e-10)
  # first SE request computes and caches the report
  s <- summary(f_lazy)
  expect_false(is.null(f_lazy$cache$sdr))
  s2 <- summary(f_eager)
  expect_equal(s$coefficients$mu, s2$coefficients$mu, tolerance = 1e-8)
})

test_that("DHARMa residuals are uniform for a correct model", {
  skip_if_not_installed("DHARMa")
  set.seed(169)
  dd <- data.frame(x = rnorm(400), g = factor(rep(1:20, 20)))
  dd$y <- rpois(400, exp(0.5 + 0.3 * dd$x + rnorm(20, 0, 0.4)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  dh <- dharma_residuals(fit, nsim = 200, seed = 1)
  expect_s3_class(dh, "DHARMa")
  ut <- DHARMa::testUniformity(dh, plot = FALSE)
  expect_gt(ut$p.value, 0.01)
  # misspecified model (ignoring the RE) should show overdispersion
  fit_bad <- frm(bf(y ~ x) + poisson(), data = dd)
  dh_bad <- dharma_residuals(fit_bad, nsim = 200, seed = 1)
  dt <- DHARMa::testDispersion(dh_bad, plot = FALSE)
  expect_lt(dt$p.value, 0.05)
})

test_that("bayesplot consumes as_tmbstan draws", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  skip_if_not_installed("bayesplot")
  set.seed(170)
  dd <- data.frame(x = rnorm(60))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  sf <- suppressWarnings(as_tmbstan(fit, chains = 1, iter = 300,
                                    refresh = 0, seed = 1))
  dr <- rstan::extract(sf, permuted = FALSE)   # iters x chains x pars
  iv <- bayesplot::mcmc_intervals_data(dr)
  expect_true(nrow(iv) >= 3)
  expect_true(all(is.finite(iv$m)))
})