# Which missing values frm() keeps latent, read from the leaves, against hand
# deletion; the window of von_mises; the warning on a read value's Laplace
# error; the methods of a fit with a reader and no latent value.

gam <- stats::Gamma(link = "log")
ll <- function(f) as.numeric(logLik(f))
# every comparison is a ratio to the fit's own size: a route equal to hand
# deletion is the same arithmetic on the same rows, up to rounding
expect_same_fit <- function(a, b) {
  expect_lt(abs(ll(a) - ll(b)) / abs(ll(b)), 1e-9)
  ea <- unlist(a$estimates[c("beta", "betad")])
  eb <- unlist(b$estimates[c("beta", "betad")])
  expect_lt(max(abs(ea - eb)) / max(abs(eb)), 1e-5)
}
pair_data <- function(seed = 12000, n = 120) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x))
  d$z <- 0.5 + 0.6 * d$y + 0.3 * d$x + rnorm(n, 0, 0.7)
  d
}

test_that("a row whose only reader has left the objective leaves it too", {
  d <- pair_data()
  m <- seq(4, nrow(d), by = 12)
  d$y[m] <- NA
  d$z[m[1:6]] <- NA
  a <- frm(bf(y | mi() ~ x) + gam + bf(z | mi() ~ mi(y) + x) + gaussian(),
           data = d)
  b <- frm(bf(y | mi() ~ x) + gam + bf(z ~ mi(y) + x) + gaussian(),
           data = d[-m[1:6], ])
  expect_setequal(a$frame$mi_map$y$unread, m[1:6])
  expect_setequal(a$frame$mi_map$z$unread, m[1:6])
  expect_length(a$estimates$miss, 4L)
  expect_same_fit(a, b)
  expect_identical(attr(logLik(a), "nobs"), attr(logLik(b), "nobs"))
})

test_that("the exgaussian the short check read as B1 again is the deletion", {
  set.seed(9503)
  n <- 200
  x <- rnorm(n)
  d <- data.frame(x = x, y = 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2))
  d$z <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
  rows <- sort(sample(n, 60))
  d$y[rows] <- NA
  d$z[rows[1:40]] <- NA
  a <- allow_warnings(frm(bf(y | mi() ~ x) + exgaussian() +
                            bf(z | mi() ~ mi(y)) + gaussian(), data = d),
                      "Laplace approximation of the")
  b <- allow_warnings(frm(bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y)) +
                            gaussian(), data = d[-rows[1:40], ]),
                      "Laplace approximation of the")
  expect_same_fit(a, b)
})

test_that("a chain of three peels from its leaf and stops at a live reader", {
  set.seed(12001)
  n <- 150
  d <- data.frame(w = rnorm(n))
  d$a <- rgamma(n, 2, 2 / exp(0.1 + 0.3 * d$w))
  d$b <- rgamma(n, 3, 3 / exp(0.1 + 0.3 * d$a))
  d$c <- 0.2 + 0.5 * d$b + rnorm(n, 0, 0.5)
  # missing at all three; at a and b; at a alone
  d$a[1:18] <- NA
  d$b[1:12] <- NA
  d$c[1:6] <- NA
  a <- allow_warnings(frm(bf(a | mi() ~ w) + gam + bf(b | mi() ~ mi(a)) + gam +
                            bf(c | mi() ~ mi(b)) + gaussian(), data = d),
                      "Laplace approximation of the")
  b <- allow_warnings(frm(bf(a | mi() ~ w) + gam + bf(b | mi() ~ mi(a)) + gam +
                            bf(c ~ mi(b)) + gaussian(), data = d[-(1:6), ]),
                      "Laplace approximation of the")
  expect_setequal(a$frame$mi_map$a$unread, 1:6)
  expect_setequal(a$frame$mi_map$b$unread, 1:6)
  expect_setequal(a$frame$mi_map$c$unread, 1:6)
  expect_length(a$frame$mi_map$a$rows, 12L)
  expect_length(a$frame$mi_map$b$rows, 6L)
  expect_same_fit(a, b)
})

test_that("a fork keeps a value while either reader stays", {
  set.seed(12002)
  n <- 150
  d <- data.frame(x = rnorm(n))
  d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x))
  d$z1 <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
  d$z2 <- -0.2 + 0.4 * d$y + rnorm(n, 0, 0.5)
  d$y[1:10] <- NA
  d$z1[1:10] <- NA
  d$z2[1:5] <- NA
  a <- frm(bf(y | mi() ~ x) + gam + bf(z1 | mi() ~ mi(y)) + gaussian() +
             bf(z2 | mi() ~ mi(y)) + gaussian(), data = d)
  b <- frm(bf(y | mi() ~ x) + gam + bf(z1 | mi() ~ mi(y)) + gaussian() +
             bf(z2 ~ mi(y)) + gaussian(), data = d[-(1:5), ])
  expect_setequal(a$frame$mi_map$y$rows, 6:10)
  expect_setequal(a$frame$mi_map$y$unread, 1:5)
  expect_length(a$frame$mi_map$z1$rows, 0L)
  expect_same_fit(a, b)
})

test_that("a reader's weight of 0 does not read", {
  d <- pair_data(12003)
  m <- seq(4, nrow(d), by = 12)
  d$y[m] <- NA
  d$wz <- ifelse(is.na(d$y), 0, 1)
  a <- frm(bf(y | mi() ~ x) + gam + bf(z | weights(wz) ~ mi(y)) + gaussian(),
           data = d)
  expect_null(a$estimates$miss)
  b <- frm(bf(y ~ x) + gam + bf(z ~ y) + gaussian(), data = d[-m, ])
  expect_lt(abs(ll(a) - ll(b)) / abs(ll(b)), 1e-9)
  d$wz[m[1:4]] <- 1
  a <- frm(bf(y | mi() ~ x) + gam + bf(z | weights(wz) ~ mi(y)) + gaussian(),
           data = d)
  b <- frm(bf(y | mi() ~ x) + gam + bf(z ~ mi(y)) + gaussian(),
           data = d[-m[5:10], ])
  expect_length(a$estimates$miss, 4L)
  expect_same_fit(a, b)
})

test_that("a reader's prediction at a row that left reads the value's mean", {
  d <- pair_data(12004)
  m <- seq(4, nrow(d), by = 12)
  d$y[m] <- NA
  d$z[m] <- NA
  fit <- frm(bf(y | mi() ~ x) + gam + bf(z | mi() ~ mi(y)) + gaussian(),
             data = d)
  fy <- fitted(fit, resp = "y")[m, "Estimate"]
  fz <- fitted(fit, resp = "z")[m, "Estimate"]
  b <- frmtmb:::fixef_by_dpar(fit)$z_mu
  eta <- b[["(Intercept)"]] + b[["miy"]] * fy
  expect_lt(max(abs(fz - eta)) / max(abs(eta)), 1e-8)
  # the value varies by row, as a placeholder would not
  expect_gt(stats::sd(fz), 0)
})

test_that("von_mises keeps its window on frm(), and a value stays inside", {
  set.seed(11400)
  n <- 150
  x <- rnorm(n)
  k <- 3
  ang <- 3.05 + 0.05 * x + stats::rnorm(n, 0, 1 / sqrt(k))
  y0 <- ((ang + pi) %% (2 * pi)) - pi
  d <- data.frame(x = x, y = y0)
  d$z <- 0.5 + 0.5 * y0 + rnorm(n, 0, 0.5)
  rows <- seq(5, n, by = 10)
  d$y[rows] <- NA
  fit <- allow_warnings(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y)) +
                              gaussian(), data = d),
                        "Laplace approximation of the")
  expect_identical(fit$frame$mi_map$y$chart, "circle")
  v <- frmtmb:::mi_values(fit, "y")[rows]
  expect_true(all(abs(v) < pi))
})

test_that("a fit whose reader skips every missing row has its methods", {
  set.seed(11600)
  n <- 120
  x <- rnorm(n)
  y0 <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * x))
  d <- data.frame(x = x, y = y0, id = seq_len(n))
  d$z <- 0.5 + 0.6 * y0 + rnorm(n, 0, 0.7)
  m <- seq(4, n, by = 12)
  d$y[m] <- NA
  d$ok <- !is.na(d$y)
  d$z2 <- ifelse(d$ok, d$z, 0)
  fit <- frm(bf(y | mi() + index(id) ~ x) + gam +
               bf(z2 | subset(ok) ~ mi(y, idx = id) + x) + gaussian(),
             data = d)
  expect_null(fit$estimates$miss)
  for (r in c("z2", "y")) {
    expect_true(all(is.finite(fitted(fit, resp = r)[, "Estimate"])), info = r)
    expect_true(all(is.finite(predict(fit, resp = r)[, "Estimate"])), info = r)
    expect_gt(NROW(residuals(fit, resp = r)), 0L)
  }
})

laplace_data <- function(slope, read, seed = 11207) {
  set.seed(seed)
  n <- 200
  x <- rnorm(n)
  y <- 0.3 + 0.5 * x - 2 + rnorm(n, 0, 0.5) + rexp(n, 1 / 2)
  d <- data.frame(x = x, y = y)
  d$z <- 0.5 + slope * y + 0.3 * x + rnorm(n, 0, 1)
  d$y[sample(n, read)] <- NA
  d
}
laplace_fit <- function(d, fam = exgaussian()) {
  seen <- FALSE
  fit <- withCallingHandlers(
    frm(bf(y | mi() ~ x) + fam + bf(z ~ mi(y) + x) + gaussian(), data = d),
    frmtmb_mi_laplace_error = function(w) {
      seen <<- TRUE
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warned = seen)
}

test_that("a weak reader of an exgaussian warns, by its own class", {
  r <- laplace_fit(laplace_data(0.05, 60))
  expect_true(r$warned)
  expect_gt(abs(r$fit$mi_laplace$tk), 1)
  expect_gte(max(r$fit$mi_laplace$moves), frmtmb:::mi_laplace_limit)
  expect_identical(r$fit$mi_laplace$n_read, 60L)
  expect_length(r$fit$mi_laplace$rows, frmtmb:::mi_laplace_max_rows)
  msg <- tryCatch(frm(bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y) + x) +
                        gaussian(), data = laplace_data(0.05, 60)),
                  frmtmb_mi_laplace_error = function(w) conditionMessage(w))
  expect_match(msg, "frm_sample()", fixed = TRUE)
  expect_match(msg, "deleting those rows", fixed = TRUE)
  expect_match(msg, "standard errors", fixed = TRUE)
  expect_match(msg, "`y_", fixed = TRUE)
})

test_that("check_mi sets what the check does", {
  d <- laplace_data(0.05, 60)
  f <- bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y) + x) + gaussian()
  off <- suppressWarnings(
    frm(f, data = d, control = frmtmb_control(check_mi = "ignore")))
  expect_null(off$mi_laplace)
  expect_error(
    suppressWarnings(
      frm(f, data = d, control = frmtmb_control(check_mi = "stop"))),
    class = "frmtmb_mi_laplace_error")
  expect_error(frmtmb_control(check_mi = "no"))
})

test_that("the parameters come back after the second stage", {
  r <- laplace_fit(laplace_data(0.05, 60))
  env <- r$fit$obj$env
  expect_equal(unname(env$last.par[env$lfixed()]),
               unname(r$fit$opt$par))
})

test_that("a von_mises value near the window's end warns by name", {
  set.seed(11600)
  n <- 200
  x <- stats::rnorm(n)
  y <- 3.14 + 0.05 * x + stats::rnorm(n, 0, 0.2)
  y <- ((y + pi) %% (2 * pi)) - pi
  d <- data.frame(x = x, y = y)
  d$z <- 0.5 + 0.7 * y + 0.3 * x + stats::rnorm(n, 0, 0.7)
  d$y[sample(n, 15)] <- NA
  said <- character()
  fit <- withCallingHandlers(
    frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(),
        data = d),
    frmtmb_mi_laplace_error = function(w) {
      said <<- c(said, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  expect_length(said, 1L)
  expect_match(said, "frm_sample()", fixed = TRUE)
  expect_match(said, "two", fixed = TRUE)
  expect_gte(fit$mi_laplace$window, 0.1)
})

test_that("a strong reader of the same family does not", {
  r <- laplace_fit(laplace_data(0.7, 60))
  expect_false(r$warned)
  ml <- r$fit$mi_laplace
  expect_true(abs(ml$tk) <= 1 ||
                max(ml$moves) < frmtmb:::mi_laplace_limit)
})

test_that("a gaussian read by a gaussian mean is exact and costs nothing", {
  set.seed(12100)
  n <- 120
  d <- data.frame(x = rnorm(n))
  d$y <- 0.3 + 0.5 * d$x + rnorm(n)
  d$z <- 0.5 + 0.7 * d$y + rnorm(n)
  d$y[1:30] <- NA
  fit <- frm(bf(y | mi() ~ x) + gaussian() + bf(z ~ mi(y)) + gaussian(),
             data = d)
  expect_null(fit$mi_laplace)
  # the same rows through the check agree: the error is rounding
  err <- frmtmb:::mi_laplace_row_errors(fit, 1:30)
  expect_lt(max(abs(err)), 1e-3)
  # a poisson reader is not a quadratic in the value
  d$k <- rpois(n, exp(0.2 + 0.2 * ifelse(is.na(d$y), 0, d$y)))
  f2 <- suppressWarnings(frm(bf(y | mi() ~ x) + gaussian() +
                               bf(k ~ mi(y)) + poisson(), data = d))
  expect_false(is.null(f2$mi_laplace))
})

test_that("importance on a fit with every missing row unread reads the rest", {
  set.seed(12005)
  n <- 240
  g <- factor(rep(1:24, each = 10))
  u <- rnorm(24, 0, 0.6)
  d <- data.frame(g = g, x = rnorm(n))
  d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x + u[g]))
  m <- seq(4, n, by = 12)
  d$y[m] <- NA
  set.seed(1)
  a <- frm(bf(y | mi() ~ x + (1 | g)), family = gam, data = d,
           importance = 200)
  set.seed(1)
  b <- frm(bf(y ~ x + (1 | g)), family = gam, data = d[-m, ], importance = 200)
  expect_lt(abs(ll(a) - ll(b)) / abs(ll(b)), 1e-3)
})

test_that("vcov(cluster =) on a fit with every missing row unread runs", {
  set.seed(12006)
  n <- 240
  g <- factor(rep(1:24, each = 10))
  d <- data.frame(g = g, x = rnorm(n))
  d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x))
  m <- seq(4, n, by = 12)
  d$y[m] <- NA
  a <- frm(bf(y | mi() ~ x), family = gam, data = d)
  b <- frm(bf(y ~ x), family = gam, data = d[-m, ])
  va <- vcov(a, cluster = d$g)
  vb <- vcov(b, cluster = d$g[-m])
  expect_lt(max(abs(va - vb)) / max(abs(vb)), 1e-6)
  # a value a density reads belongs to no cluster
  d$z <- rnorm(n)
  r <- suppressWarnings(frm(bf(y | mi() ~ x) + gam + bf(z ~ mi(y)) + gaussian(),
                            data = d))
  expect_error(vcov(r, cluster = d$g), "latent values")
})

test_that("mi_sampling() names what it takes", {
  d <- pair_data(12007)
  d$y[1:5] <- NA
  fit <- frm(bf(y | mi() ~ x), family = gam, data = d)
  expect_error(mi_sampling(fit, "bad"), "must be one of")
})

test_that("the laplace seam gives an unfitted object the fit's objective", {
  d <- pair_data(12008)
  m <- seq(4, nrow(d), by = 12)
  d$y[m] <- NA
  o <- frm(bf(y | mi() ~ x), family = gam, data = d, dry_run = "objective")
  a <- mi_sampling(o, "laplace")
  expect_null(a$estimates$miss)
  expect_setequal(a$frame$mi_map$y$unread, m)
  b <- frm(bf(y ~ x), family = gam, data = d[-m, ], dry_run = "objective")
  expect_lt(abs(a$obj$fn(b$obj$par) - b$obj$fn(b$obj$par)) /
              abs(b$obj$fn(b$obj$par)), 1e-10)
})

test_that("df.residual() and anova() count the rows a missing y left", {
  d <- pair_data(12009)
  m <- seq(4, nrow(d), by = 12)
  d$y[m] <- NA
  a <- frm(bf(y | mi() ~ x), family = gam, data = d)
  b <- frm(bf(y ~ x), family = gam, data = d[-m, ])
  b0 <- frm(bf(y ~ 1), family = gam, data = d[-m, ])
  # the mi() fit's likelihood is on the complete rows, so every count
  # agrees with the complete-case fit's: through 0.81.0 df.residual()
  # counted the rows mi() left unread, and anova() refused to compare
  # the mi() fit with a complete-case fit
  expect_identical(attr(logLik(a), "nobs"), nrow(d) - length(m))
  expect_equal(df.residual(a), df.residual(b))
  ref <- anova(b, b0)$Chisq[2L]
  expect_no_error(tab <- anova(a, b0))
  expect_lt(abs(tab$Chisq[2L] - ref), 1e-6 * abs(ref))
})
