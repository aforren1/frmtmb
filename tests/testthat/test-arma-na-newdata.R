# predict(newdata = ) under cov = FALSE with the response NA in some
# rows. brms's .predictor_arma() runs the recursion in each group's
# time order and fills a missing response with a draw from the family
# at the row's shifted mean, so the rows after it read that draw's
# residual (brms 2.23.0 source, dev/formrobust-findings.md section 3).
# Before, frmtmb refused the call by name (ledger row
# brmsfit-methods:747).

an_fit <- local({
  set.seed(31)
  G <- 30
  Tn <- 8
  d <- expand.grid(t = 1:Tn, g = factor(1:G))
  d$x <- stats::rnorm(nrow(d))
  e <- as.vector(apply(matrix(stats::rnorm(G * Tn), Tn, G), 2, function(z) {
    as.vector(stats::filter(z, 0.6, "recursive"))
  }))
  d$y <- 1 + 0.5 * d$x + e
  list(d = d, fit = frm(bf(y ~ x + ar(t, g, p = 1)), data = d))
})

an_newdata <- function() {
  nd <- an_fit$d[an_fit$d$g %in% c("1", "2"), ]
  nd$y[nd$g == "1" & nd$t >= 5] <- NA
  nd
}

test_that("fitted() fills a missing response with its expected value", {
  fit <- an_fit$fit
  nd <- an_newdata()
  f <- fitted(fit, newdata = nd)
  expect_false(anyNA(f[, "Estimate"]))
  ar1 <- unname(fit$estimates[["thetaac"]][1L])
  b <- fixef(fit)[, "Estimate"]
  mu <- b[["Intercept"]] + b[["x"]] * nd$x
  i4 <- which(nd$g == "1" & nd$t == 4)
  # row 5 reads the observed residual of row 4; row 6 reads row 5's
  # fill, which is its one-step mean, so its residual is the shift
  m5 <- mu[i4 + 1L] + ar1 * (nd$y[i4] - mu[i4])
  m6 <- mu[i4 + 2L] + ar1 * (m5 - mu[i4 + 1L])
  expect_lt(abs(f[i4 + 1L, "Estimate"] - m5) / abs(m5),
            64 * .Machine$double.eps)
  expect_lt(abs(f[i4 + 2L, "Estimate"] - m6) / abs(m6),
            64 * .Machine$double.eps)
  # rows with every earlier response observed are what they were
  full <- an_fit$d[an_fit$d$g %in% c("1", "2"), ]
  ff <- fitted(fit, newdata = full)
  ok <- !(nd$g == "1" & nd$t >= 6)
  expect_identical(f[ok, "Estimate"], ff[ok, "Estimate"])
})

test_that("predict() fills a missing response with a draw, as brms", {
  fit <- an_fit$fit
  nd <- an_newdata()
  set.seed(1)
  p <- predict(fit, newdata = nd, ndraws = 4000, propagate_error = FALSE)
  expect_false(anyNA(p[, "Estimate"]))
  ar1 <- unname(fit$estimates[["thetaac"]][1L])
  sig <- unname(sigma(fit)[1L])
  i6 <- which(nd$g == "1" & nd$t == 6)
  i5 <- i6 - 1L
  # row 5 is drawn around a mean every input of which is observed; row
  # 6 adds the variance of row 5's fill through the ar coefficient, so
  # its sd is sigma sqrt(1 + ar^2). The expected-value fill would give
  # sigma, a ratio of 1 / sqrt(1 + ar^2) = 0.88 here.
  r5 <- p[i5, "Est.Error"] / sig
  r6 <- p[i6, "Est.Error"] / (sig * sqrt(1 + ar1^2))
  # 4000 draws: the sd of a sample sd ratio is about 1 / sqrt(8000)
  expect_lt(abs(r5 - 1), 4 / sqrt(8000))
  expect_lt(abs(r6 - 1), 4 / sqrt(8000))
  # the guard-absent case: every response observed draws row 6 around
  # the one-step mean with sd sigma
  full <- an_fit$d[an_fit$d$g %in% c("1", "2"), ]
  set.seed(1)
  pf <- predict(fit, newdata = full, ndraws = 4000, propagate_error = FALSE)
  expect_lt(abs(pf[i6, "Est.Error"] / sig - 1), 4 / sqrt(8000))
})

test_that("newdata without the response column is a forecast", {
  fit <- an_fit$fit
  nd <- an_newdata()
  nd$y <- NULL
  set.seed(2)
  p <- predict(fit, newdata = nd, ndraws = 50)
  expect_identical(dim(p), c(nrow(nd), 4L))
  expect_false(anyNA(p))
  expect_false(anyNA(fitted(fit, newdata = nd)))
  # a response of the wrong type is still refused by name
  nd$y <- "a"
  expect_error(predict(fit, newdata = nd, ndraws = 5),
               "The response must be numeric", fixed = TRUE)
})

test_that("an arma() fit with an NA response in newdata predicts", {
  set.seed(9)
  d <- an_fit$d
  fit <- frm(bf(y ~ x + arma(t, g, p = 1, q = 1)), data = d)
  nd <- d[1:10, ]
  nd$y[8:10] <- NA
  set.seed(3)
  p <- predict(fit, newdata = nd, ndraws = 20)
  expect_false(anyNA(p[, "Estimate"]))
})

test_that("an arma() fill leaves the residual net of its MA term", {
  # brms's .predictor_arma() stores err = y - eta_before_ar, where
  # eta_before_ar already holds the MA term; with the fill at its
  # expected value the residual of a filled row is then its AR term
  # alone. A residual taken as y - mu, the MA term kept in, moves every
  # later row by ma * (that MA term).
  d <- an_fit$d
  fit <- frm(bf(y ~ x + arma(t, g, p = 1, q = 1)), data = d)
  th <- unname(fit$estimates[["thetaac"]])
  ar1 <- th[1L]
  ma1 <- th[2L]
  b <- fixef(fit)[, "Estimate"]
  nd <- d[d$g == "1", ]
  nd$y[5:6] <- NA
  mu <- b[["Intercept"]] + b[["x"]] * nd$x
  m <- numeric(nrow(nd))
  err <- 0
  for (t in seq_len(nrow(nd))) {
    sma <- if (t > 1L) ma1 * err else 0
    sar <- if (t > 1L) ar1 * err else 0
    m[t] <- mu[t] + sma + sar
    yt <- if (is.na(nd$y[t])) m[t] else nd$y[t]
    err <- yt - mu[t] - sma
  }
  f <- fitted(fit, newdata = nd)[, "Estimate"]
  expect_lt(max(abs(f - m) / abs(m)), 64 * .Machine$double.eps)
})
