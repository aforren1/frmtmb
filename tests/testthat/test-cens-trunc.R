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
  expect_vector_equal(fixef_by_dpar(fit)$mu, unname(coef(ref)), tol = 1e-4)
  expect_lt(abs(exp(fixef_by_dpar(fit)$sigma[[1]]) - ref$scale), 1e-3)
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
  expect_vector_equal(fixef_by_dpar(fit)$mu, c(1, 0.5), tol = 0.15)
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
    expect_vector_equal(fixef_by_dpar(fit)$mu, fixef_by_dpar(fit_num)$mu,
                        tol = 1e-8)
    expect_vector_equal(fixef_by_dpar(fit)$sigma, fixef_by_dpar(fit_num)$sigma,
                        tol = 1e-8)
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
  expect_vector_equal(fixef_by_dpar(fit_s)$mu, fixef_by_dpar(fit_n)$mu,
                      tol = 1e-8)
})

test_that("cens/trunc validation", {
  dd <- data.frame(y = rpois(50, 3), x = rnorm(50), cen = 0)
  # a discrete family is admitted on the CDF it supplies, so poisson
  # and negbinomial fit and nbinom1, which carries none, is refused for
  # that (brms has no nbinom1 to say otherwise)
  expect_s3_class(frm(bf(y | cens(cen) ~ x) + poisson(), data = dd),
                  "frmtmb_fit")
  expect_s3_class(frm(bf(y | cens(cen) ~ x) + negbinomial(), data = dd),
                  "frmtmb_fit")
  expect_error(frm(bf(y | cens(cen) ~ x) + nbinom1(), data = dd),
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
# Discrete censoring and truncation are brms 2.23.0's, term for term:
# right censoring at k is Y > k, left censoring at k is Y <= k, an
# interval is k < Y <= k2, and only trunc(lb = ) is inclusive, through
# F(lb - 1). Through frmtmb 0.72.0 every censoring lower edge was
# inclusive instead.
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

# brms's reading; `incl = TRUE` is the inclusive one this replaced
disc_cens_ll <- function(d, m, incl = FALSE) {
  s <- if (incl) 1 else 0
  ll <- numeric(nrow(d))
  i_o <- d$cc == "none"
  i_r <- d$cc == "right"
  i_l <- d$cc == "left"
  i_i <- d$cc == "interval"
  ll[i_o] <- stats::dpois(d$y[i_o], m, log = TRUE)
  ll[i_r] <- stats::ppois(d$y[i_r] - s, m, lower.tail = FALSE,
                          log.p = TRUE)
  ll[i_l] <- stats::ppois(d$y[i_l], m, log.p = TRUE)
  ll[i_i] <- log(stats::ppois(d$y2[i_i], m) -
                   stats::ppois(d$y[i_i] - s, m))
  sum(ll)
}

test_that("all four discrete censoring codes match a hand-rolled likelihood", {
  d <- disc_cens_data()
  fit <- frm(y | cens(cc, y2) ~ 1, family = poisson(), data = d)
  mu <- exp(unname(fixef_by_dpar(fit)$mu))
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
    mk <- exp(unname(fixef_by_dpar(f)$mu))
    expect_equal(as.numeric(logLik(f)), disc_cens_ll(dk, mk),
                 tolerance = 1e-10, label = k)
  }
})

test_that("the discrete convention is brms's, not the inclusive one", {
  # the two readings differ by the point mass at the recorded value,
  # which is large enough to move the estimate
  d <- disc_cens_data()
  d <- d[d$cc %in% c("none", "right"), ]
  fit <- frm(y | cens(cc) ~ 1, family = poisson(), data = d)
  mu <- exp(unname(fixef_by_dpar(fit)$mu))
  expect_equal(as.numeric(logLik(fit)), disc_cens_ll(d, mu),
               tolerance = 1e-10)
  expect_gt(disc_cens_ll(d, mu, incl = TRUE) - as.numeric(logLik(fit)), 1)
})

test_that("a right-censored count under trunc(lb = ) alone is brms's row", {
  d <- disc_cens_data()
  d <- d[d$cc %in% c("none", "right") & d$y >= 2, ]
  fit <- frm(y | cens(cc) + trunc(lb = 2) ~ 1, family = poisson(),
             data = d)
  mu <- exp(unname(fixef_by_dpar(fit)$mu))
  i_r <- d$cc == "right"
  ll <- numeric(nrow(d))
  ll[!i_r] <- stats::dpois(d$y[!i_r], mu, log = TRUE)
  ll[i_r] <- stats::ppois(d$y[i_r], mu, lower.tail = FALSE, log.p = TRUE)
  # the window is P(Y >= 2), F(lb - 1) for the inclusive lower bound
  ref <- sum(ll) - nrow(d) * stats::ppois(1, mu, lower.tail = FALSE,
                                          log.p = TRUE)
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-10)
})

test_that("a one-point discrete interval is an event of probability zero", {
  set.seed(708)
  d <- data.frame(y = stats::rpois(200, 4))
  d$y2 <- d$y
  d$cc <- "interval"
  expect_error(frm(y | cens(cc, y2) ~ 1, family = poisson(), data = d),
               "must exceed the lower bounds")
  dc <- data.frame(y = rnorm(50), cc = "interval")
  dc$y2 <- dc$y
  expect_error(frm(y | cens(cc, y2) ~ 1, family = gaussian(), data = dc),
               "must exceed the lower bounds")
})

test_that("a fractional discrete bound is refused by name", {
  set.seed(709)
  d <- data.frame(y = stats::rpois(60, 4), cc = "none", y2 = NA_real_)
  d$cc[1:10] <- "interval"
  d$y2[1:10] <- d$y[1:10] + 2.5
  expect_error(frm(y | cens(cc, y2) ~ 1, family = poisson(), data = d),
               "must be a whole number")
  d <- data.frame(y = stats::rpois(60, 4) + 2)
  expect_error(frm(y | trunc(lb = 1.5) ~ 1, family = poisson(), data = d),
               "must be a whole number")
})

test_that("a custom discrete family with a CDF is censored too", {
  # what admits a family is the CDF it supplies, not the type it
  # declares; this one has a probability-scale lcdf alone
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
  d$cc <- ifelse(d$y > 7, "right", "none")
  d$y <- pmin(d$y, 7)
  fit <- frm(y | cens(cc) ~ 1, family = poisson(), data = d)
  expect_error(residuals(fit, type = "osa"),
               "not supported on a cens.. fit with a discrete family")
  # every other residual type still works
  expect_true(all(is.finite(residuals(fit, type = "response")[, "Estimate"])))
})

test_that("a count left-censored at the support minimum is P(Y = 0)", {
  # left censoring at 0 is P(Y <= 0), the mass at zero itself; brms
  # reads it the same way
  set.seed(711)
  d <- data.frame(y = stats::rpois(120, 2))
  d$cc <- ifelse(d$y == 0, "left", "none")
  expect_gt(sum(d$cc == "left"), 5)
  fit <- frm(y | cens(cc) ~ 1, family = poisson(), data = d)
  plain <- frm(y ~ 1, family = poisson(), data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(plain)),
               tolerance = 1e-10)
})

test_that("trunc(lb <= 0) on a count is no truncation, as in brms", {
  set.seed(712)
  d <- data.frame(x = stats::rnorm(150))
  d$y <- stats::rpois(150, exp(0.8 + 0.3 * d$x))
  plain <- frm(y ~ x, family = poisson(), data = d)
  zero <- frm(y | trunc(lb = 0) ~ x, family = poisson(), data = d)
  expect_equal(as.numeric(logLik(zero)), as.numeric(logLik(plain)),
               tolerance = 1e-10)
  # per row: 0 and -2 are absent, 1 is P(Y >= 1), and an infinite lower
  # bound is absent too
  d$lb <- ifelse(d$y >= 1, rep_len(c(-2, 0, 1), 150), 0)
  d$ub <- Inf
  fit <- frm(y | trunc(lb = lb, ub = ub) ~ x, family = poisson(), data = d)
  mu <- frmtmb:::eval_dpars(fit)[[1]]$mu
  ref <- sum(stats::dpois(d$y, mu, log = TRUE) -
               ifelse(d$lb >= 1, stats::ppois(d$lb - 1, mu,
                                              lower.tail = FALSE,
                                              log.p = TRUE), 0))
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-10)
  inf <- frm(y | trunc(lb = -Inf, ub = 30) ~ x, family = poisson(),
             data = d)
  mu <- frmtmb:::eval_dpars(inf)[[1]]$mu
  expect_equal(as.numeric(logLik(inf)),
               sum(stats::dpois(d$y, mu, log = TRUE) -
                     stats::ppois(30, mu, log.p = TRUE)),
               tolerance = 1e-10)
})

test_that("an infinite truncation bound is absent and fits", {
  # through frmtmb 0.72.0 trunc(lb = -Inf, ub = Inf)
  # stopped on a NaN gradient from dnorm(-Inf) * Inf
  set.seed(713)
  d <- data.frame(x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(120)
  plain <- frm(y ~ x, family = gaussian(), data = d)
  open <- frm(y | trunc(lb = -Inf, ub = Inf) ~ x, family = gaussian(),
              data = d)
  expect_equal(as.numeric(logLik(open)), as.numeric(logLik(plain)),
               tolerance = 1e-10)
})

test_that("the trunc() window check refuses lb >= ub and admits the bounds", {
  set.seed(714)
  d <- data.frame(y = stats::rpois(80, 4) + 1)
  expect_error(frm(y | trunc(lb = 5, ub = 5) ~ 1, family = poisson(),
                   data = d[d$y == 5, , drop = FALSE]),
               "lb >= ub")
  expect_error(frm(y | trunc(lb = 6, ub = 2) ~ 1, family = gaussian(),
                   data = data.frame(y = c(3, 4, 5))),
               "lb >= ub")
  # a response ON either bound is inside the window, on a count and on a
  # continuous response, and the check is absent there
  lo <- min(d$y)
  hi <- max(d$y)
  fit <- frm(y | trunc(lb = lo, ub = hi) ~ 1, family = poisson(),
             data = cbind(d, lo = lo, hi = hi))
  expect_s3_class(fit, "frmtmb_fit")
  mu <- exp(unname(fixef_by_dpar(fit)$mu))
  expect_equal(as.numeric(logLik(fit)),
               sum(stats::dpois(d$y, mu, log = TRUE)) -
                 nrow(d) * log(stats::ppois(hi, mu) - stats::ppois(lo - 1, mu)),
               tolerance = 1e-10)
  # enough rows that the window's two parameters are well determined:
  # on five rows the fit stopped at false convergence under OpenBLAS
  dg <- data.frame(y = c(-1, 2, seq(-0.9, 1.9, length.out = 40)))
  fg <- frm(y | trunc(lb = -1, ub = 2) ~ 1, family = gaussian(), data = dg)
  expect_equal(fg$opt$convergence, 0L)
  expect_error(frm(y | trunc(lb = 0) ~ 1, family = gaussian(), data = dg),
               paste0("outside of the truncation bounds: ", sum(dg$y < 0),
                      " of ", nrow(dg), " rows"))
})

test_that("discrete censoring gives no notice: it is brms's reading", {
  # the once-per-session notice of the inclusive reading is gone, with
  # the divergence it announced
  withr::local_options(frmtmb.notices = TRUE)
  local_frmtmb_registries()
  st <- frmtmb:::frmtmb_notice_state
  rm(list = ls(envir = st, all.names = TRUE), envir = st)
  d <- disc_cens_data()
  got <- character(0)
  withCallingHandlers(
    frm(y | cens(cc, y2) ~ 1, family = poisson(), data = d),
    message = function(m) {
      got <<- c(got, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_length(got, 0)
})

test_that("a censored row under both trunc() bounds is the event inside the window", {
  # frmtmb's composition, not brms 2.23.0's: brms divides P(Y > y) or
  # P(Y <= y), the event over the whole line, by the window's mass
  # (brms-34), and differs by 6.16 and 14.85 on a poisson and a
  # negative binomial in trunc(lb = 2, ub = 9)
  set.seed(716)
  d <- data.frame(y = stats::rpois(300, 4))
  d <- d[d$y >= 2 & d$y <= 9, , drop = FALSE]
  # a right-censored row at ub would be the empty event (9, 9]
  d$cc <- ifelse(d$y %in% 7:8, 1, ifelse(d$y <= 3, -1, 0))
  fit <- frm(y | cens(cc) + trunc(lb = 2, ub = 9) ~ 1, family = poisson(),
             data = d)
  mu <- exp(unname(fixef_by_dpar(fit)$mu))
  Fp <- function(q) stats::ppois(q, mu)
  Z <- Fp(9) - Fp(1)
  ll <- ifelse(d$cc == 1, log(Fp(9) - Fp(d$y)),
               ifelse(d$cc == -1, log(Fp(d$y) - Fp(1)),
                      stats::dpois(d$y, mu, log = TRUE))) - log(Z)
  expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-10)
  whole <- ifelse(d$cc == 1, stats::ppois(d$y, mu, lower.tail = FALSE,
                                          log.p = TRUE),
                  ifelse(d$cc == -1, log(Fp(d$y)),
                         stats::dpois(d$y, mu, log = TRUE))) - log(Z)
  expect_gt(sum(whole) - sum(ll), 1)
  # and on a continuous response
  set.seed(717)
  g <- data.frame(y = stats::rnorm(400, 1, 1.5))
  g <- g[g$y > -1 & g$y < 3, , drop = FALSE]
  g$cc <- ifelse(g$y > 2.5, 1, ifelse(g$y < -0.5, -1, 0))
  g$y <- pmin(pmax(g$y, -0.5), 2.5)
  fg <- frm(y | cens(cc) + trunc(lb = -1, ub = 3) ~ 1, family = gaussian(),
            data = g)
  m <- unname(fixef_by_dpar(fg)$mu)
  s <- exp(unname(fixef_by_dpar(fg)$sigma))
  P <- function(q) stats::pnorm(q, m, s)
  Zg <- P(3) - P(-1)
  lg <- ifelse(g$cc == 1, log(P(3) - P(g$y)),
               ifelse(g$cc == -1, log(P(g$y) - P(-1)),
                      stats::dnorm(g$y, m, s, log = TRUE))) - log(Zg)
  expect_equal(as.numeric(logLik(fg)), sum(lg), tolerance = 1e-10)
})

test_that("a count censored right at its trials is refused by name", {
  # P(Y > n) = 0 for a count out of n trials; brms 2.23.0 accepts the
  # row, its log density is -Inf, and it makes no draw
  set.seed(718)
  d <- data.frame(y = stats::rbinom(40, 5, 0.4), tr = 5, cc = 0)
  d$y[1] <- 5
  d$cc[1] <- 1
  expect_error(frm(y | trials(tr) + cens(cc) ~ 1, family = binomial(),
                   data = d), "probability zero")
  expect_error(frm(y | trials(tr) + cens(cc) ~ 1, family = beta_binomial(),
                   data = d), "probability zero")
  # left censoring at 0 is P(Y = 0), which is possible, and fits
  d$y[1] <- 0
  d$cc[1] <- -1
  fit <- frm(y | trials(tr) + cens(cc) ~ 1, family = binomial(), data = d)
  expect_true(is.finite(as.numeric(logLik(fit))))
})
