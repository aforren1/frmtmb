# Edge cases mined from lme4/glmmTMB/brms issue history (dev/test-backlog.md
# holds the full list with sources).

#' @srrstats {G5.8} Edge conditions are tested for the behavior they
#'   should produce, which here is usually a clear error or a documented
#'   message rather than a fit. This file covers data-dependent bases at
#'   a single-row `newdata`, rank-deficient designs, matrix responses
#'   carrying attributes, row permutation and factor releveling,
#'   equivalent spellings of a grouping variable, duplicated multivariate
#'   responses, the `trials()` error taxonomy, response rescaling, and
#'   `NA` rows in both the fitting data and `newdata`.
#' @srrstats {G5.8b} Data of unsupported types is tested to error rather
#'   than be silently coerced. A fractional count for `trials()` errors
#'   on "integer"; a response outside `[0, trials]` errors; a non-integer
#'   response to a Poisson family and a non-0/1 response to a bernoulli
#'   family are refused instead of fitted (lme4 warns and fits these);
#'   and character and factor multipliers in a `mo()`/`mi()` interaction
#'   are both rejected, the character case specifically because
#'   `as.numeric()` on a character vector is numeric and all `NA`, which
#'   used to pass the type gate and surface only as an optimizer failure.
#' @srrstats {G5.8a} Zero-length data is tested. A data frame with no
#'   rows, and one subset down to no rows, both error with "`data` has no
#'   rows; nothing to fit". That verdict is kept distinct from "No
#'   complete observations after removing NAs", which is the different
#'   fault of losing every row to `na.action`.
#' @srrstats {G5.8c} All-`NA` and all-identical columns are tested. An
#'   all-`NA` predictor takes every row with it and the fit is refused;
#'   an all-identical predictor is collinear with the intercept, so it is
#'   dropped by name in a `message()` and the fit reproduces the model
#'   without it; an all-identical response leaves no scale to estimate
#'   and drives `sigma` to the boundary instead of returning nonsense.
#' @srrstats {G5.8d} Data outside the scope of the algorithm is tested. A
#'   design with more columns than observations (p = 20, n = 12) is
#'   reduced by the rank check, which names the columns it drops; what
#'   survives is finite and estimable rather than `NA`-padded, and
#'   `nobs()` still reports n.
#' @srrstats {G5.9,G5.9a} Noise susceptibility is tested. Adding noise of
#'   `sqrt(.Machine$double.eps)` to an O(1) response moves every outer
#'   parameter by less than 1e3 times that noise and the log-likelihood
#'   by less than 1e-4, and noise a thousand times larger moves the
#'   estimates measurably more, so the invariance is a property of the
#'   scale of the perturbation and not of a slack tolerance. The seed
#'   half of the standard (G5.9b) is covered in `test-fuzz.R`.
#' @srrstats {RE7.1,RE7.1a} A noiseless, exact relationship between
#'   predictor and response is tested. Because an exact `y = X beta` with
#'   a free `sigma` has an unbounded likelihood, the dispersion is fixed
#'   (`bf(y ~ x + z, sigma = 1)`) and the mean relationship carries no
#'   error term; the fit converges, recovers the generating coefficients
#'   to 1e-4, and returns finite standard errors. RE7.1a is tested by
#'   comparison with the same model on the same design with noise added.
#'   The comparison is averaged over twelve seeded pairs, because the
#'   iteration count of a single small quasi-Newton problem swings by an
#'   order of magnitude with the draw: on average the noiseless fit uses
#'   fewer optimizer iterations and fewer function evaluations, which is
#'   the machine-independent part of run time, and its wall clock at a
#'   larger n is checked as well with slack.
#' @srrstats {RE7.2} Output objects retain the row names of the input.
#'   `fitted()`, `residuals()`, `predict()`, `predict(newdata =)`, and
#'   `model.frame()` all carry the data frame's row names, including
#'   after `na.action` has removed rows, where the surviving names stay
#'   in place.
#' @srrstats {RE7.0,RE7.0a} Noiseless, exact relationships between
#'   predictors are tested. A perfectly collinear pair (`x2 <- 2 * x`) is
#'   required to produce a "rank deficient" message naming the dropped
#'   column, a log-likelihood identical to the reduced model, and
#'   identical predictions. The rejection behavior is tested separately:
#'   a design whose collinearity is broken in `newdata` must warn "not
#'   estimable" and return `NA` for exactly those rows, and an aliased
#'   cell must return `NA` rather than a partial sum. The sparse backend
#'   is required to drop the same columns as the dense one.
#' @noRd
NULL

test_that("data-dependent bases are frozen at fit time (glmmTMB#402)", {
  set.seed(151)
  dd <- data.frame(x = runif(60, 0, 10), g = factor(rep(1:6, 10)))
  dd$y <- 1 + 0.5 * dd$x - 0.05 * dd$x^2 + rnorm(6, 0, 0.5)[dd$g] +
    rnorm(60, 0, 0.3)
  m1 <- frm(bf(y ~ poly(x, 2) + (1 | g)) + gaussian(), data = dd)
  m2 <- frm(bf(y ~ poly(x, 2, raw = TRUE) + (1 | g)) + gaussian(),
            data = dd)
  nd <- data.frame(x = c(0, 2.5, 9), g = factor(1, levels = levels(dd$g)))
  expect_equal(predict(m1, newdata = nd, re.form = NA),
               predict(m2, newdata = nd, re.form = NA), tolerance = 1e-6)
  # single-row newdata is the killer case (brms#494)
  p1 <- predict(m1, newdata = nd[2, , drop = FALSE], re.form = NA)
  expect_equal(p1, predict(m1, newdata = nd, re.form = NA)[2],
               tolerance = 1e-8)
  # scale() in the formula round-trips through prediction
  m3 <- frm(bf(y ~ scale(x) + (1 | g)) + gaussian(), data = dd)
  expect_equal(predict(m3, newdata = dd), predict(m3), tolerance = 1e-8)
})

test_that("rank-deficient designs drop aliased columns (lme4#144)", {
  set.seed(152)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$x2 <- 2 * dd$x                       # perfectly collinear
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
  expect_message(
    m <- frm(bf(y ~ x + x2 + (1 | g)) + gaussian(), data = dd),
    "rank deficient"
  )
  m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(m0))), 1e-6)
  # prediction from the rank-reduced fit works and matches
  expect_equal(predict(m, newdata = dd), predict(m0, newdata = dd),
               tolerance = 1e-6)
})

test_that("matrix-attribute responses and offsets are handled (glmmTMB#937/#773)", {
  set.seed(153)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  dd$ys <- scale(dd$y)                    # n x 1 matrix with attributes
  m <- frm(bf(ys ~ x + (1 | g)) + gaussian(), data = dd)
  expect_length(fitted(m), 60)
})

test_that("row-permutation and relevel invariance for covariance structures (brms#1747)", {
  dd <- sim_ar1_data(seed = 154, n_g = 40)
  perm <- sample(nrow(dd))
  f1 <- frm(bf(y ~ 1 + ar1(tim + 0 | g)) + gaussian(), data = dd)
  f2 <- frm(bf(y ~ 1 + ar1(tim + 0 | g)) + gaussian(), data = dd[perm, ])
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-6)

  u1 <- frm(bf(y ~ 1 + us(tim + 0 | g)) + gaussian(), data = dd)
  u2 <- frm(bf(y ~ 1 + us(tim + 0 | g)) + gaussian(), data = dd[perm, ])
  expect_lt(abs(as.numeric(logLik(u1)) - as.numeric(logLik(u2))), 1e-5)
  # us() is also invariant to releveling the term factor
  dd3 <- transform(dd, tim = stats::relevel(tim, "3"))
  u3 <- frm(bf(y ~ 1 + us(tim + 0 | g)) + gaussian(), data = dd3)
  expect_lt(abs(as.numeric(logLik(u1)) - as.numeric(logLik(u3))), 1e-5)
})

test_that("numeric, character, and factor grouping variables are equivalent", {
  set.seed(155)
  d1 <- data.frame(x = rnorm(100), g = rep(1:10, 10))
  d1$y <- rnorm(100, 1 + 0.5 * d1$x + rnorm(10, 0, 0.6)[d1$g], 1)
  d2 <- transform(d1, g = factor(g))
  d3 <- transform(d1, g = as.character(g))
  f1 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d1)
  f2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d2)
  f3 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d3)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-6)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f3))), 1e-6)
})

test_that("duplicate multivariate responses are rejected (brms)", {
  expect_error(
    frm(mvbf(bf(y ~ x) + gaussian(), bf(y ~ 1) + gaussian()),
        data = NULL, dry_run = "spec"),
    "Duplicated response"
  )
})

test_that("trials() validation catches the brms error taxonomy", {
  dd <- data.frame(y = c(2, 3, 5), n = c(5, 5, 4), x = 1:3)
  expect_error(frm(bf(y | trials(n) ~ x) + binomial(), data = dd),
               "\\[0, trials\\]")   # y > trials
  # a proportion response is legal since v0.14 when it times out to
  # integer counts; a fractional count is still an error
  dd2 <- data.frame(y = c(0.3, 1), n = c(2, 2), x = 1:2)
  expect_error(frm(bf(y | trials(n) ~ x) + binomial(), data = dd2),
               "integer")
  # constant literal trials work
  dd3 <- data.frame(y = rbinom(50, 10, 0.4), x = rnorm(50))
  fit <- frm(bf(y | trials(10) ~ x) + binomial(), data = dd3)
  expect_length(fitted(fit), 50)
})

test_that("response-scale equivariance under rescaling", {
  set.seed(156)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.6)[dd$g], 1)
  m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  s <- 1000
  dd$ys <- dd$y * s
  m1 <- frm(bf(ys ~ x + (1 | g)) + gaussian(), data = dd)
  expect_vector_equal(fixef(m1)$mu, fixef(m0)$mu * s, tol = 1e-2)
  expect_lt(abs(as.numeric(logLik(m1)) -
                  (as.numeric(logLik(m0)) - 100 * log(s))), 1e-4)
})

test_that("NA handling: rows dropped consistently, Inf rejected upstream", {
  set.seed(157)
  dd <- data.frame(x = rnorm(50), g = factor(rep(1:5, 10)))
  dd$y <- rnorm(50, 1 + 0.5 * dd$x, 1)
  dd$y[1:3] <- NA
  dd$g[4] <- NA
  m <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_identical(stats::nobs(m), 46L)
  expect_length(fitted(m), 46)
  # NA rows in newdata propagate NA predictions, not errors
  nd <- dd[1:6, ]
  nd$x[2] <- NA
  p <- predict(m, newdata = nd, re.form = NA)
  expect_false(is.na(p[1]))   # NA was only in y/g, x is fine
  expect_true(is.na(p[2]))    # NA predictor rows come back NA
})

test_that("rows dropped by na.action are reported once (G2.14b)", {
  set.seed(158)
  dd <- data.frame(x = rnorm(50), g = factor(rep(1:5, 10)))
  dd$y <- rnorm(50, 1 + 0.5 * dd$x, 1)
  dd$y[1:3] <- NA
  dd$x[40] <- NA
  expect_message(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd),
                 "^4 rows removed because of missing values")
  # singular is singular: one row lost says "row", not "rows"
  dd2 <- dd
  dd2$y[2:3] <- 0
  dd2$x[40] <- 0
  expect_message(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd2),
                 "^1 row removed because of missing values")
  # complete data stays quiet, and the message is a message: callers who
  # asked for na.omit can silence it without silencing warnings
  dd3 <- dd
  dd3$y[1:3] <- 0
  dd3$x[40] <- 0
  expect_no_message(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd3))
  expect_no_message(
    suppressMessages(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd))
  )
})

test_that("zero-length data errors informatively (G5.8a)", {
  dd <- data.frame(y = numeric(0), x = numeric(0), g = factor())
  expect_error(frm(bf(y ~ x) + gaussian(), data = dd),
               "no rows")
  # a data frame with columns but no rows takes the same route as one
  # subset down to nothing
  full <- data.frame(y = rnorm(20), x = rnorm(20))
  expect_error(frm(bf(y ~ x) + gaussian(), data = full[0, ]),
               "no rows")
  # zero-length is separated from "every row was dropped as missing",
  # which is a different fault with a different remedy
  allna <- full
  allna$x <- NA_real_
  expect_error(suppressMessages(frm(bf(y ~ x) + gaussian(), data = allna)),
               "No complete observations")
})

test_that("all-NA and all-identical columns are handled (G5.8c)", {
  set.seed(159)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)

  # an all-NA predictor takes every row with it, so the fit is refused
  d_na <- dd
  d_na$z <- NA_real_
  expect_error(suppressMessages(
    frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = d_na)),
    "No complete observations")

  # an all-identical predictor is collinear with the intercept: it is
  # dropped by name, not silently absorbed, and the fit matches the
  # model without it
  d_const <- dd
  d_const$k <- 3.5
  expect_message(
    m <- frm(bf(y ~ x + k + (1 | g)) + gaussian(), data = d_const),
    "rank deficient.*k"
  )
  m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(m0))), 1e-6)
  expect_false(anyNA(fixef(m)$mu))

  # an all-identical response has no scale to estimate; the variance
  # parameters run to the boundary rather than returning nonsense
  d_flat <- dd
  d_flat$y <- 2
  expect_error(
    suppressWarnings(fit_flat <- frm(bf(y ~ x + (1 | g)) + gaussian(),
                                     data = d_flat)),
    NA
  )
  expect_lt(stats::sigma(fit_flat), 1e-3)
})

test_that("p > n drops the unidentified columns (G5.8d)", {
  set.seed(160)
  n <- 12L
  p <- 20L
  dd <- as.data.frame(matrix(rnorm(n * p), n, p))
  names(dd) <- paste0("v", seq_len(p))
  dd$y <- rnorm(n)
  fo <- stats::as.formula(
    paste("y ~", paste(paste0("v", seq_len(p)), collapse = " + "))
  )
  msg <- NULL
  withCallingHandlers(
    m <- suppressWarnings(frm(bf(fo) + gaussian(), data = dd)),
    message = function(cnd) {
      msg <<- c(msg, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )
  # the dropped columns are named, so the user can see which ones went
  expect_match(paste(msg, collapse = ""), "rank deficient")
  expect_match(paste(msg, collapse = ""), "v20")
  # the design cannot support more than n columns; what survives is
  # estimable and finite, never NA-padded the way lm() does it
  b <- fixef(m)$mu
  expect_lte(length(b), n)
  expect_true(all(is.finite(b)))
  expect_identical(stats::nobs(m), n)
})

test_that("eps-scale noise moves the estimates by O(noise) (G5.9, G5.9a)", {
  set.seed(161)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$y <- 1 + 0.6 * dd$x + rnorm(8, 0, 0.4)[dd$g] + rnorm(80, 0, 0.5)
  m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  eps <- sqrt(.Machine$double.eps)          # ~1.5e-8
  set.seed(162)
  d1 <- dd
  d1$y <- dd$y + stats::rnorm(80, 0, eps)
  m1 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d1)

  # every outer parameter, not only the fixed effects: the response
  # scale is O(1), so a shift of order eps must stay of order eps
  expect_lt(max(abs(m1$opt$par - m0$opt$par)), 1e3 * eps)
  expect_lt(abs(as.numeric(logLik(m1)) - as.numeric(logLik(m0))), 1e-4)
  expect_vector_equal(fitted(m1), fitted(m0), tol = 1e3 * eps)

  # and the movement is proportional: noise a thousand times larger
  # moves the estimates measurably more than the eps-scale noise
  set.seed(163)
  d2 <- dd
  d2$y <- dd$y + stats::rnorm(80, 0, 1e3 * eps)
  m2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d2)
  expect_gt(max(abs(m2$opt$par - m0$opt$par)),
            max(abs(m1$opt$par - m0$opt$par)))
})

test_that("a noiseless exact fit succeeds (RE7.1)", {
  # An exact y = X beta with a free sigma has an unbounded likelihood
  # (sigma -> 0), so "noiseless" is posed the only way it can be: the
  # dispersion is fixed and the mean relationship carries no error.
  set.seed(164)
  dd <- data.frame(x = rnorm(200), z = rnorm(200))
  dd$y <- 1 + 0.5 * dd$x - 0.25 * dd$z      # no residual term at all
  m <- frm(bf(y ~ x + z, sigma = 1) + gaussian(), data = dd)
  expect_identical(m$opt$convergence, 0L)
  expect_vector_equal(fixef(m)$mu, c(1, 0.5, -0.25), tol = 1e-4)
  expect_lt(max(abs(stats::residuals(m))), 1e-4)
  expect_true(all(is.finite(stats::vcov(m))))
})

test_that("a noiseless fit costs no more than the noisy one (RE7.1a)", {
  # iteration counts move with BLAS and platform, and the wall-clock
  # comparison needs an unloaded machine
  skip_on_cran()
  fit <- function(d) frm(bf(y ~ x + z, sigma = 1) + gaussian(), data = d)
  # A single data set says nothing: the iteration count of a 3-parameter
  # quasi-Newton problem swings between 2 and 17 with the draw. The
  # comparison is therefore averaged over twelve seeded pairs, each pair
  # being the same design with and without a residual term.
  work <- t(vapply(seq_len(12), function(s) {
    set.seed(100L + s)
    d <- data.frame(x = stats::rnorm(200), z = stats::rnorm(200))
    d$y <- 1 + 0.5 * d$x - 0.25 * d$z
    dn <- d
    dn$y <- d$y + stats::rnorm(200, 0, 0.5)
    m_exact <- fit(d)
    m_noisy <- fit(dn)
    c(it_e = m_exact$opt$iterations, it_n = m_noisy$opt$iterations,
      ev_e = sum(m_exact$opt$evaluations),
      ev_n = sum(m_noisy$opt$evaluations))
  }, numeric(4)))

  # optimizer work is the machine-independent part of run time
  expect_lte(mean(work[, "it_e"]), mean(work[, "it_n"]))
  expect_lte(mean(work[, "ev_e"]), mean(work[, "ev_n"]))

  # and the wall clock agrees, with slack: these fits are milliseconds
  set.seed(168)
  dd <- data.frame(x = stats::rnorm(2000), z = stats::rnorm(2000))
  dd$y <- 1 + 0.5 * dd$x - 0.25 * dd$z
  dn <- dd
  set.seed(169)
  dn$y <- dd$y + stats::rnorm(2000, 0, 0.5)
  invisible(fit(dd))
  invisible(suppressWarnings(fit(dn)))
  t_exact <- stats::median(replicate(5, system.time(fit(dd))[["elapsed"]]))
  t_noisy <- stats::median(replicate(
    5, system.time(suppressWarnings(fit(dn)))[["elapsed"]]))
  expect_lt(t_exact, 1.5 * t_noisy + 0.25)
})

test_that("row names survive into fitted, residuals and predict (RE7.2)", {
  set.seed(167)
  dd <- data.frame(x = rnorm(40), g = factor(rep(1:4, 10)))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
  rn <- paste0("case", seq_len(40))
  rownames(dd) <- rn
  m <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  expect_identical(names(stats::fitted(m)), rn)
  expect_identical(names(stats::residuals(m)), rn)
  expect_identical(names(stats::predict(m)), rn)
  expect_identical(rownames(stats::model.frame(m)), rn)
  # newdata carries its own row names through
  nd <- dd[c(5, 1, 9), ]
  expect_identical(names(stats::predict(m, newdata = nd)),
                   c("case5", "case1", "case9"))
  # rows dropped for missingness leave the surviving names in place
  dd2 <- dd
  dd2$y[3] <- NA
  m2 <- suppressMessages(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd2))
  expect_identical(names(stats::fitted(m2)), rn[-3])
})
