test_that("right-censored gaussian matches survival::survreg", {
  skip_if_not_installed("survival")
  set.seed(101)
  n <- 400
  x <- rnorm(n)
  ystar <- 1 + 0.8 * x + rnorm(n, 0, 1.2)
  cpoint <- 2
  y <- pmin(ystar, cpoint)
  cen <- as.numeric(ystar > cpoint)   # 1 = right-censored
  dd <- data.frame(y = y, x = x, cen = cen)

  fit <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd)
  ref <- survival::survreg(survival::Surv(y, 1 - cen) ~ x,
                           data = dd, dist = "gaussian")
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-5)
  expect_vector_equal(fixef(fit)$mu, unname(coef(ref)), tol = 1e-4)
  expect_lt(abs(exp(fixef(fit)$sigma[[1]]) - ref$scale), 1e-3)
})

test_that("mixed left/right censoring matches a hand-rolled reference", {
  set.seed(102)
  n <- 400
  x <- rnorm(n)
  ystar <- 0.5 * x + rnorm(n)
  lo <- -1; hi <- 1.5
  y <- pmin(pmax(ystar, lo), hi)
  cen <- ifelse(ystar < lo, -1, ifelse(ystar > hi, 1, 0))
  dd <- data.frame(y = y, x = x, cen = cen)
  fit <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd)

  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * x
    s <- exp(p$ls)
    ll <- RTMB::dnorm(y, mu, s, log = TRUE)
    Fv <- RTMB::pnorm((y - mu) / s)
    ir <- which(cen == 1); il <- which(cen == -1)
    ll[ir] <- log(1 - Fv[ir])
    ll[il] <- log(Fv[il])
    -sum(ll)
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("truncation matches a hand-rolled reference", {
  set.seed(103)
  n <- 2000
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n)
  keep <- y > 0
  dd <- data.frame(y = y[keep], x = x[keep])

  fit <- frm(bf(y | trunc(lb = 0) ~ x) + gaussian(), data = dd)
  yv <- dd$y; xv <- dd$x
  nll_ref <- function(p) {
    mu <- p$b[1] + p$b[2] * xv
    s <- exp(p$ls)
    -sum(RTMB::dnorm(yv, mu, s, log = TRUE) -
           log(1 - RTMB::pnorm((0 - mu) / s)))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  # truncation-corrected fit recovers the latent coefficients
  expect_vector_equal(fixef(fit)$mu, c(1, 0.5), tol = 0.15)
})

test_that("brms character and factor censoring codes fit identically", {
  set.seed(104)
  n <- 300
  x <- rnorm(n)
  ystar <- 0.5 * x + rnorm(n)
  lo <- -1; hi <- 1.5
  y <- pmin(pmax(ystar, lo), hi)
  cen <- ifelse(ystar < lo, -1, ifelse(ystar > hi, 1, 0))
  lab <- c("left", "none", "right")[cen + 2]
  dd <- data.frame(y = y, x = x, cen = cen, lab = lab,
                   flab = factor(lab),
                   # prefix matching, case-insensitively
                   pre = c("LEFT", "n", "r")[cen + 2])

  fit_num <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd)
  for (nm in c("lab", "flab", "pre")) {
    fo <- stats::as.formula(paste0("y | cens(", nm, ") ~ x"))
    fit <- frm(bf(fo) + gaussian(), data = dd)
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(fit_num))),
              1e-8)
    expect_vector_equal(fixef(fit)$mu, fixef(fit_num)$mu, tol = 1e-8)
    expect_vector_equal(fixef(fit)$sigma, fixef(fit_num)$sigma, tol = 1e-8)
  }
})

test_that("string interval censoring matches its numeric code", {
  set.seed(105)
  n <- 200
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n)
  iv <- y > 0.5
  dd <- data.frame(y = y, x = x, y2 = y + 1,
                   cs = ifelse(iv, "interval", "none"),
                   cn = ifelse(iv, 2, 0))
  fit_s <- frm(bf(y | cens(cs, y2) ~ x) + gaussian(), data = dd)
  fit_n <- frm(bf(y | cens(cn, y2) ~ x) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit_s)) - as.numeric(logLik(fit_n))), 1e-8)
  expect_vector_equal(fixef(fit_s)$mu, fixef(fit_n)$mu, tol = 1e-8)
})

test_that("cens/trunc validation", {
  dd <- data.frame(y = rpois(50, 3), x = rnorm(50), cen = 0)
  # a discrete family is admitted on the CDF it supplies, so poisson
  # fits and negbinomial, which carries none, is refused for that
  expect_s3_class(frm(bf(y | cens(cen) ~ x) + poisson(), data = dd),
                  "frmtmb_fit")
  expect_error(frm(bf(y | cens(cen) ~ x) + negbinomial(), data = dd),
               "family with a CDF")
  dd2 <- data.frame(y = rnorm(50), x = rnorm(50), cen = 3)
  expect_error(frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd2),
               "codes")
  dd3 <- data.frame(y = rnorm(50), x = rnorm(50), cen = 2)
  expect_error(frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd3),
               "needs upper bounds")
  expect_error(frm(bf(y | trunc(0) ~ x) + gaussian(),
                   data = NULL, dry_run = "spec"),
               "named bounds")
})

test_that("unusable censoring labels error informatively", {
  dd <- data.frame(y = rnorm(50), x = rnorm(50),
                   bad = rep(c("right", "sideways"), 25),
                   amb = "",
                   iv = "interval")
  expect_error(frm(bf(y | cens(bad) ~ x) + gaussian(), data = dd),
               "sideways")
  expect_error(frm(bf(y | cens(bad) ~ x) + gaussian(), data = dd),
               "interval")
  # an empty label prefix-matches every spelling, so it is not a code
  expect_error(frm(bf(y | cens(amb) ~ x) + gaussian(), data = dd),
               "cannot decode")
  # "interval" still needs the upper bounds, same as the numeric code
  expect_error(frm(bf(y | cens(iv) ~ x) + gaussian(), data = dd),
               "needs upper bounds")
  # factor labels go through the same decoder
  dd$fbad <- factor(dd$bad)
  expect_error(frm(bf(y | cens(fbad) ~ x) + gaussian(), data = dd),
               "sideways")
})

# ---------------------------------------------------------------------
# Discrete censoring. A bound NAMES a value the response can take and is
# INCLUDED: right censoring at k is Y >= k, an interval is
# k <= Y <= k2, and left censoring at k is Y <= k. Equivalently every
# LOWER edge enters the CDF as F(edge - 1), which is the shift
# trunc(lb = ) has always applied.
# ---------------------------------------------------------------------

disc_cens_data <- function(seed = 707, n = 240, lam = 5) {
  set.seed(seed)
  ytrue <- stats::rpois(n, lam)
  code <- rep(c("none", "right", "left", "interval"), length.out = n)
  y <- ytrue
  y2 <- rep(NA_real_, n)
  i_r <- code == "right"
  i_l <- code == "left"
  i_i <- code == "interval"
  y[i_r] <- pmin(ytrue[i_r], 7)
  y[i_l] <- pmax(ytrue[i_l], 3)
  y[i_i] <- pmax(ytrue[i_i] - 1, 0)
  y2[i_i] <- y[i_i] + 2
  data.frame(y = y, y2 = y2, cc = code)
}

disc_cens_ll <- function(d, m) {
  ll <- numeric(nrow(d))
  i_o <- d$cc == "none"
  i_r <- d$cc == "right"
  i_l <- d$cc == "left"
  i_i <- d$cc == "interval"
  ll[i_o] <- stats::dpois(d$y[i_o], m, log = TRUE)
  ll[i_r] <- stats::ppois(d$y[i_r] - 1, m, lower.tail = FALSE,
                          log.p = TRUE)
  ll[i_l] <- stats::ppois(d$y[i_l], m, log.p = TRUE)
  ll[i_i] <- log(stats::ppois(d$y2[i_i], m) -
                   stats::ppois(d$y[i_i] - 1, m))
  sum(ll)
}

test_that("all four discrete censoring codes match a hand-rolled likelihood", {
  d <- disc_cens_data()
  fit <- frm(y | cens(cc, y2) ~ 1, family = poisson(), data = d)
  mu <- exp(unname(fixef(fit)$mu))
  expect_equal(as.numeric(logLik(fit)), disc_cens_ll(d, mu),
               tolerance = 1e-10)
  # and the estimate is the hand-rolled likelihood's own maximum
  o <- stats::optimize(function(g) -disc_cens_ll(d, exp(g)), c(-2, 4),
                       tol = 1e-12)
  expect_equal(mu, exp(o$minimum), tolerance = 1e-6)

  # each code alone, against the same reference
  for (k in c("right", "left", "interval")) {
    dk <- d[d$cc %in% c("none", k), ]
    f <- frm(y | cens(cc, y2) ~ 1, family = poisson(), data = dk)
    mk <- exp(unname(fixef(f)$mu))
    expect_equal(as.numeric(logLik(f)), disc_cens_ll(dk, mk),
                 tolerance = 1e-10, label = k)
  }
})

test_that("the discrete convention is NOT the continuous one", {
  # the difference is the point mass at the recorded value, and it is
  # large enough to move the estimate: this is what the old blanket
  # refusal of discrete families was standing in for
  d <- disc_cens_data()
  d <- d[d$cc %in% c("none", "right"), ]
  fit <- frm(y | cens(cc) ~ 1, family = poisson(), data = d)
  mu <- exp(unname(fixef(fit)$mu))
  i_r <- d$cc == "right"
  excl <- sum(stats::dpois(d$y[!i_r], mu, log = TRUE)) +
    sum(stats::ppois(d$y[i_r], mu, lower.tail = FALSE, log.p = TRUE))
  expect_equal(as.numeric(logLik(fit)), disc_cens_ll(d, mu),
               tolerance = 1e-10)
  expect_gt(as.numeric(logLik(fit)) - excl, 1)
})

test_that("discrete cens() composes with trunc() on one convention", {
  d <- disc_cens_data()
  d <- d[d$cc %in% c("none", "right") & d$y >= 2, ]
  fit <- frm(y | cens(cc) + trunc(lb = 2) ~ 1, family = poisson(),
             data = d)
  mu <- exp(unname(fixef(fit)$mu))
  i_r <- d$cc == "right"
  ll <- numeric(nrow(d))
  ll[!i_r] <- stats::dpois(d$y[!i_r], mu, log = TRUE)
  ll[i_r] <- stats::ppois(d$y[i_r] - 1, mu, lower.tail = FALSE,
                          log.p = TRUE)
  # one window, one rule: F(lb - 1) for the normalizer as for the rows
  ref <- sum(ll) - nrow(d) * log(1 - stats::ppois(1, mu))
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-10)
})

test_that("a one-point discrete interval is the exact observation", {
  set.seed(708)
  d <- data.frame(y = stats::rpois(200, 4))
  d$y2 <- d$y
  d$cc <- "interval"
  band <- frm(y | cens(cc, y2) ~ 1, family = poisson(), data = d)
  plain <- frm(y ~ 1, family = poisson(), data = d)
  expect_equal(as.numeric(logLik(band)), as.numeric(logLik(plain)),
               tolerance = 1e-10)
  # on a continuous response the same interval has probability zero
  dc <- data.frame(y = rnorm(50), cc = "interval")
  dc$y2 <- dc$y
  expect_error(frm(y | cens(cc, y2) ~ 1, family = gaussian(), data = dc),
               "must exceed the lower bounds")
})

test_that("a non-integer discrete censoring bound is refused by name", {
  set.seed(709)
  d <- data.frame(y = stats::rpois(60, 4), cc = "none")
  d$cc[1:10] <- "right"
  d$y <- d$y + 0.5
  expect_error(frm(y | cens(cc) ~ 1, family = poisson(), data = d),
               "support is the integers")
})

test_that("a custom discrete family with a CDF is censored too", {
  # the point of the change: what admits a family is the CDF it
  # supplies, not the type it declares
  pois_cdf <- frmtmb_family(
    "pois_cdf", dpars = "mu", links = list(mu = "log"),
    type = "discrete", accepts_aterms = c("weights", "cens", "trunc"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, dpars[["mu"]], log = TRUE)
    },
    lcdf = function(q, dpars, aterms) RTMB::ppois(q, dpars[["mu"]]),
    init_dpars = list(mu = function(y, aterms) mean(y)))
  no_cdf <- frmtmb_family(
    "pois_no_cdf", dpars = "mu", links = list(mu = "log"),
    type = "discrete",
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, dpars[["mu"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y)))

  d <- disc_cens_data()
  fit <- frm(bf(y | cens(cc, y2) ~ 1) + pois_cdf, data = d)
  ref <- frm(y | cens(cc, y2) ~ 1, family = poisson(), data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-10)
  expect_error(frm(bf(y | cens(cc, y2) ~ 1) + no_cdf, data = d),
               "family with a CDF")
})

test_that("osa residuals are refused on a censored discrete fit", {
  set.seed(710)
  d <- data.frame(y = stats::rpois(200, 4))
  d$cc <- ifelse(d$y >= 7, "right", "none")
  d$y <- pmin(d$y, 7)
  fit <- frm(y | cens(cc) ~ 1, family = poisson(), data = d)
  expect_error(residuals(fit, type = "osa"),
               "not supported on a cens.. fit with a discrete family")
  # every other residual type still works
  expect_true(all(is.finite(residuals(fit, type = "response"))))
})

test_that("a discrete row censored at the support minimum is free", {
  # right censoring at 0 says "0 or more", which is P(Y >= 0) = 1 and
  # no information at all. The shift reads F(-1), so this is also the
  # test that the CDF is evaluated below the support without harm.
  set.seed(711)
  d <- data.frame(y = stats::rpois(120, 2))
  d$cc <- ifelse(d$y == 0, "right", "none")
  expect_gt(sum(d$cc == "right"), 5)
  fit <- suppressWarnings(frm(y | cens(cc) ~ 1, family = poisson(),
                              data = d))
  mu <- exp(unname(fixef(fit)$mu))
  obs <- d$cc == "none"
  expect_equal(as.numeric(logLik(fit)),
               sum(stats::dpois(d$y[obs], mu, log = TRUE)),
               tolerance = 1e-10)
})
