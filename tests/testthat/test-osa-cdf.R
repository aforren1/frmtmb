# One-step-ahead residuals read off the family's CDF (osa_method = "cdf",
# the default for a core family with a CDF, but an untruncated gaussian()).
#
# Reference: with no group-level term the one-step predictive of a row is
# its fitted distribution, so the residual of a count lies in
# [qnorm(F(y - 1)), qnorm(F(y))] for every valid randomization, and the
# residual of a continuous response is qnorm(F(y)). F is formed here from
# stats:: functions at the fitted parameters.

osa_r <- function(fit, ...) residuals(fit, type = "osa", ...)[, "Estimate"]
osa_dp <- function(fit) {
  lapply(frmtmb:::eval_dpars(fit)[[1L]], as.numeric)
}

# TRUE where r lies in [qnorm(lo), qnorm(hi)], with a slack of 1e-8 of the
# row's own interval: every valid randomization passes and nothing else
osa_in_step <- function(r, lo, hi) {
  qlo <- stats::qnorm(lo)
  qhi <- stats::qnorm(pmin(hi, 1))
  wid <- qhi - qlo
  wid[!is.finite(wid)] <- 1
  is.finite(r) & r >= qlo - 1e-8 * wid & r <= qhi + 1e-8 * wid
}

osa_design <- function(seed = 9, n = 150) {
  set.seed(seed)
  x <- stats::rnorm(n)
  p <- stats::plogis(-0.3 + 0.8 * x)
  list(x = x, p = p, n = n)
}

test_that("count residuals are randomized quantile residuals of the fit", {
  # PIN: on 0.79.1 oneStepGeneric with no discreteSupport integrated a
  # spline through the mass function's continuous extension. Here
  # 3 poisson, 3 binomial, 53 bernoulli and 8 negbinomial rows were NaN
  # or outside their own step, and beta_binomial was
  # refused (dev/osa-findings.md). Seen failing there.
  z <- osa_design()
  k <- 10L
  pbb <- function(q, a, b) {
    vapply(seq_along(q), function(i) {
      if (q[i] < 0) return(0)
      j <- 0:q[i]
      sum(exp(lchoose(k, j) + lbeta(j + a[i], k - j + b[i]) -
                lbeta(a[i], b[i])))
    }, 0)
  }
  mc <- exp(0.8 + 0.5 * z$x)
  cases <- list(
    poisson = list(y = stats::rpois(z$n, mc), fam = poisson(),
                   F = function(q, dp) stats::ppois(q, dp$mu)),
    binomial = list(y = stats::rbinom(z$n, k, z$p), fam = binomial(),
                    F = function(q, dp) stats::pbinom(q, k, dp$mu)),
    bernoulli = list(y = stats::rbinom(z$n, 1, z$p), fam = bernoulli(),
                     F = function(q, dp) stats::pbinom(q, 1, dp$mu)),
    negbinomial = list(y = stats::rnbinom(z$n, size = 2, mu = mc),
                       fam = negbinomial(),
                       F = function(q, dp) {
                         stats::pnbinom(q, size = dp$shape, mu = dp$mu)
                       }),
    beta_binomial = list(y = stats::rbinom(z$n, k, stats::rbeta(
      z$n, z$p * 8, (1 - z$p) * 8)), fam = beta_binomial(),
      F = function(q, dp) {
        pbb(q, dp$mu * dp$phi, (1 - dp$mu) * dp$phi)
      }))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    d <- data.frame(y = cs$y, x = z$x, k = k)
    form <- if (nm %in% c("binomial", "beta_binomial")) {
      bf(y | trials(k) ~ x)
    } else bf(y ~ x)
    fit <- frm(form, family = cs$fam, data = d)
    r <- osa_r(fit)
    dp <- osa_dp(fit)
    ok <- osa_in_step(r, cs$F(d$y - 1, dp), cs$F(d$y, dp))
    expect_identical(sum(!ok), 0L, info = nm)
  }
})

test_that("continuous residuals are qnorm of the fitted CDF", {
  # PIN: on 0.79.1 Gamma read 0.043 and Beta 0.105 off qnorm(F(y)),
  # relative to the largest residual, and lognormal stopped with an S4
  # type error; on dev/osa-audit.R's design inverse.gaussian read NaN on
  # 3 rows with no warning. Seen failing there.
  z <- osa_design()
  mp <- exp(0.3 + 0.4 * z$x)
  cases <- list(
    Gamma = list(y = stats::rgamma(z$n, 3, 3 / mp), fam = Gamma(link = "log"),
                 F = function(y, dp) {
                   stats::pgamma(y, dp$shape, dp$shape / dp$mu)
                 }),
    beta = list(y = stats::rbeta(z$n, z$p * 6, (1 - z$p) * 6),
                fam = Beta(),
                F = function(y, dp) {
                  stats::pbeta(y, dp$mu * dp$phi, (1 - dp$mu) * dp$phi)
                }),
    lognormal = list(y = exp(0.2 + 0.4 * z$x + stats::rnorm(z$n, 0, 0.5)),
                     fam = lognormal(),
                     F = function(y, dp) stats::plnorm(y, dp$mu, dp$sigma)),
    student = list(y = 1 + 0.7 * z$x + stats::rt(z$n, 5), fam = student(),
                   F = function(y, dp) {
                     stats::pt((y - dp$mu) / dp$sigma, dp$nu)
                   }),
    inverse.gaussian = list(
      y = stats::rgamma(z$n, 3, 3 / mp), fam = inverse.gaussian(link = "log"),
      F = function(y, dp) {
        s <- sqrt(dp$shape / y)
        stats::pnorm(s * (y / dp$mu - 1)) +
          exp(2 * dp$shape / dp$mu) * stats::pnorm(-s * (y / dp$mu + 1))
      }))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    d <- data.frame(y = cs$y, x = z$x)
    fit <- frm(bf(y ~ x), family = cs$fam, data = d)
    r <- osa_r(fit)
    ref <- stats::qnorm(cs$F(d$y, osa_dp(fit)))
    expect_true(all(is.finite(r)), info = nm)
    expect_lt(max(abs(r - ref)) / max(abs(ref)), 1e-8)
  }
})

test_that("a zero-inflated or hurdle count reads its point mass as data", {
  # PIN: refused by name on 0.79.1, because the density branches on
  # y == 0 and the generic route hands it no value to compare
  z <- osa_design()
  mc <- exp(0.8 + 0.5 * z$x)
  zi <- stats::runif(z$n) < 0.25
  d <- data.frame(x = z$x, yz = ifelse(zi, 0, stats::rpois(z$n, mc)),
                  yh = ifelse(zi, 0, 1 + stats::rpois(z$n, mc)))
  fz <- frm(bf(yz ~ x), family = zero_inflated_poisson(), data = d)
  dp <- osa_dp(fz)
  Fz <- function(q) ifelse(q < 0, 0, dp$zi + (1 - dp$zi) *
                             stats::ppois(q, dp$mu))
  expect_identical(sum(!osa_in_step(osa_r(fz), Fz(d$yz - 1), Fz(d$yz))), 0L)
  fh <- frm(bf(yh ~ x), family = hurdle_poisson(), data = d)
  dp <- osa_dp(fh)
  Fh <- function(q) {
    pos <- (stats::ppois(q, dp$mu) - stats::dpois(0, dp$mu)) /
      (1 - stats::dpois(0, dp$mu))
    ifelse(q < 0, 0, dp$hu + (1 - dp$hu) * pmax(pos, 0))
  }
  expect_identical(sum(!osa_in_step(osa_r(fh), Fh(d$yh - 1), Fh(d$yh))), 0L)
  # the other routes still hand the density the osa object, and say so
  expect_error(osa_r(fz, osa_method = "oneStepGeneric", discreteSupport = 0:30),
               "osa_method = \"cdf\" reads the count as data", fixed = TRUE)
})

test_that("a truncated count is read inside its window, row by row", {
  # PIN: on 0.79.1 trunc(lb = 1) put 4 of 150 rows outside their step,
  # and a row-varying bound was refused (dev/osa-edge.R)
  z <- osa_design()
  lb <- rep(c(1, 2), length.out = z$n)
  m0 <- exp(1.2 + 0.4 * z$x)
  y <- stats::qpois(stats::runif(z$n, stats::ppois(lb - 1, m0), 1), m0)
  d <- data.frame(x = z$x, y = y, lb = lb)
  for (vary in c(TRUE, FALSE)) {
    dv <- d
    if (!vary) {
      dv$lb <- 1
      dv$y <- pmax(dv$y, 1)
    }
    fit <- frm(y | trunc(lb = lb) ~ x, family = poisson(), data = dv)
    m <- osa_dp(fit)$mu
    Z <- 1 - stats::ppois(dv$lb - 1, m)
    lo <- (stats::ppois(dv$y - 1, m) - stats::ppois(dv$lb - 1, m)) / Z
    hi <- (stats::ppois(dv$y, m) - stats::ppois(dv$lb - 1, m)) / Z
    expect_identical(sum(!osa_in_step(osa_r(fit), lo, hi)), 0L,
                     info = paste("row-varying:", vary))
  }
})

test_that("a truncated gaussian with a group-level term has every row", {
  # PIN: on 0.79.1 oneStepGeneric integrated each row with the later
  # rows' normalizers left in. Here the integral failed on every row, and
  # oneStepPredict() stopped with "$ operator is invalid for atomic
  # vectors" (dev/osa-trunc-re.R)
  set.seed(31)
  d <- data.frame(g = factor(rep(1:12, each = 8)), x = stats::rnorm(96))
  u <- stats::rnorm(12, 0, 0.8)[d$g]
  d$y <- NA_real_
  todo <- seq_len(96)
  while (length(todo)) {
    d$y[todo] <- 1 + 0.5 * d$x[todo] + u[todo] + stats::rnorm(length(todo))
    todo <- which(d$y <= 0.5)
  }
  fit <- frm(y | trunc(lb = 0.5) ~ x + (1 | g), family = gaussian(),
             data = d)
  r <- osa_r(fit)
  expect_true(all(is.finite(r)))
  expect_lt(abs(mean(r)) / (stats::sd(r) / sqrt(length(r))), 4)
  # against the exact one-step predictive at the fit's estimates: the
  # truncated-normal CDF averaged over the group effect given the earlier
  # rows, on a grid. TMB's Laplace error was 1.6e-3 here
  # (dev/osa-trunc-re.R); a route that ignores the window is far off
  dp <- osa_dp(fit)
  eta <- dp$mu - fit$estimates$b[as.integer(d$g)]
  sig <- dp$sigma[1]
  sdu <- exp(fit$estimates$theta[[1]])
  ug <- seq(-8, 8, length.out = 801) * sdu
  ex <- numeric(96)
  for (gi in split(seq_len(96), d$g)) {
    w <- stats::dnorm(ug, 0, sdu)
    for (i in gi) {
      m <- eta[i] + ug
      Z <- stats::pnorm((m - 0.5) / sig)
      Fy <- (stats::pnorm((d$y[i] - m) / sig) -
               stats::pnorm((0.5 - m) / sig)) / Z
      ex[i] <- stats::qnorm(sum(w * Fy) / sum(w))
      w <- w * stats::dnorm(d$y[i], m, sig) / Z
    }
  }
  expect_lt(max(abs(r - ex)) / max(abs(ex)), 0.01)
})

# Where a residual sits in its own step, recovered from two calls with
# different seeds: pnorm(r) = hi - U p for each call's uniforms U. The
# upper end hi of the step is the one-step CDF at y, whatever the seed.
osa_step_hi <- function(fit) {
  ra <- osa_r(fit, seed = 123)
  rb <- osa_r(fit, seed = 456)
  n <- length(ra)
  ua <- local({ set.seed(123); stats::runif(n) })
  ub <- local({ set.seed(456); stats::runif(n) })
  p <- (stats::pnorm(ra) - stats::pnorm(rb)) / (ub - ua)
  stats::pnorm(ra) + ua * p
}

osa_re_counts <- function() {
  set.seed(301)
  d <- data.frame(g = factor(rep(1:20, each = 10)), t = rep(1:10, 20),
                  x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(0.2 + 0.4 * d$x +
                                 stats::rnorm(20, 0, 0.7)[d$g]))
  d
}

test_that("a weight is a frequency under a group-level term", {
  # PIN: with the weights left out of the one-step tape, the later rows
  # were conditioned on an unweighted posterior of the group effect, and
  # a weight of 2 differed from the row written twice by up to 0.33
  # (dev/reviews/2026-11-06-osa.md, B1). Weight 2 on rows t <= 3 is the
  # same fitted model as those rows duplicated, so every later row has
  # the same one-step CDF in both
  d <- osa_re_counts()
  d$w <- ifelse(d$t <= 3, 2, 1)
  fw <- frm(y | weights(w) ~ x + (1 | g), family = poisson(), data = d)
  dd <- rbind(d, d[d$t <= 3, ])
  dd <- dd[order(dd$g, dd$t), ]
  fd <- frm(y ~ x + (1 | g), family = poisson(), data = dd)
  hw <- osa_step_hi(fw)[d$t >= 4]
  hd <- osa_step_hi(fd)[dd$t >= 4]
  expect_lt(max(abs(stats::qnorm(hw) - stats::qnorm(hd))) /
              max(abs(stats::qnorm(hd))), 1e-4)
})

test_that("a row of weight 0 is NA and leaves the other rows alone", {
  # PIN: weight 0 took the rows out of the fit but not out of the
  # one-step conditioning, and the mean residual of the other rows read
  # -1.751 against 0.014 without them (dev/reviews/2026-11-06-osa.md)
  d <- osa_re_counts()
  d$y[d$t <= 2] <- d$y[d$t <= 2] + 6
  d$w <- ifelse(d$t <= 2, 0, 1)
  fw <- frm(y | weights(w) ~ x + (1 | g), family = poisson(), data = d)
  fd <- frm(y ~ x + (1 | g), family = poisson(), data = d[d$t > 2, ])
  expect_true(all(is.na(osa_r(fw)[d$t <= 2])))
  hw <- osa_step_hi(fw)[d$t > 2]
  hd <- osa_step_hi(fd)
  expect_lt(max(abs(stats::qnorm(hw) - stats::qnorm(hd))) /
              max(abs(stats::qnorm(hd))), 1e-4)
})

test_that("an outlier's residual is read in its own tail", {
  # PIN: TMB forms qnorm(F - U p) on the probability scale, so an upper
  # tail below about 1e-16 read Inf, with a warning that the CDF could
  # not be formed (dev/reviews/2026-11-06-osa.md, M4)
  z <- osa_design()
  d <- data.frame(x = z$x, y = stats::rpois(z$n, exp(1 + 0.2 * z$x)))
  d$y[30] <- 60
  fit <- frm(y ~ x, family = poisson(), data = d)
  r <- osa_r(fit)
  m <- osa_dp(fit)$mu
  # the step on the upper side, from the log upper tails at y - 1 and y
  qlo <- -stats::qnorm(stats::ppois(d$y - 1, m, lower.tail = FALSE,
                                    log.p = TRUE), log.p = TRUE)
  qhi <- -stats::qnorm(stats::ppois(d$y, m, lower.tail = FALSE,
                                    log.p = TRUE), log.p = TRUE)
  expect_true(is.finite(r[30]))
  expect_gt(r[30], qlo[30] - 1e-8 * (qhi[30] - qlo[30]))
  expect_lt(r[30], qhi[30] + 1e-8 * (qhi[30] - qlo[30]))
  dg <- data.frame(x = z$x, y = 1 + 0.5 * z$x + stats::rnorm(z$n))
  dg$y[7] <- dg$y[7] + 12
  dg$y[8] <- dg$y[8] - 12
  fg <- frm(y ~ x, family = gaussian(), data = dg)
  rg <- osa_r(fg, osa_method = "cdf")
  dpg <- osa_dp(fg)
  zg <- (dg$y - dpg$mu) / dpg$sigma
  expect_lt(max(abs(rg - zg)) / max(abs(zg)), 1e-10)
})

test_that("the zero-inflated and hurdle counts sit in their own steps", {
  z <- osa_design()
  k <- 10L
  mc <- exp(0.8 + 0.5 * z$x)
  zi <- stats::runif(z$n) < 0.25
  d <- data.frame(
    x = z$x, k = k,
    ynb = ifelse(zi, 0, stats::rnbinom(z$n, size = 2, mu = mc)),
    yhn = ifelse(zi, 0, 1 + stats::rnbinom(z$n, size = 2, mu = mc)),
    yb = ifelse(zi, 0, stats::rbinom(z$n, k, z$p)),
    ybb = ifelse(zi, 0, stats::rbinom(z$n, k, stats::rbeta(
      z$n, z$p * 8, (1 - z$p) * 8))))
  pbb <- function(q, a, b) {
    vapply(seq_along(q), function(i) {
      if (q[i] < 0) return(0)
      j <- 0:q[i]
      sum(exp(lchoose(k, j) + lbeta(j + a[i], k - j + b[i]) -
                lbeta(a[i], b[i])))
    }, 0)
  }
  zif <- function(q, zi, Fb) ifelse(q < 0, 0, zi + (1 - zi) * Fb)
  cases <- list(
    list(bf(ynb ~ x), zero_inflated_negbinomial(), "ynb", function(q, dp) {
      zif(q, dp$zi, stats::pnbinom(q, size = dp$shape, mu = dp$mu))
    }),
    list(bf(yhn ~ x), hurdle_negbinomial(), "yhn", function(q, dp) {
      p0 <- stats::dnbinom(0, size = dp$shape, mu = dp$mu)
      pos <- (stats::pnbinom(q, size = dp$shape, mu = dp$mu) - p0) / (1 - p0)
      ifelse(q < 0, 0, dp$hu + (1 - dp$hu) * pmax(pos, 0))
    }),
    list(bf(yb | trials(k) ~ x), zero_inflated_binomial(), "yb",
         function(q, dp) zif(q, dp$zi, stats::pbinom(q, k, dp$mu))),
    list(bf(ybb | trials(k) ~ x), zero_inflated_beta_binomial(), "ybb",
         function(q, dp) {
           zif(q, dp$zi, pbb(q, dp$mu * dp$phi, (1 - dp$mu) * dp$phi))
         }))
  for (cs in cases) {
    fit <- allow_warnings(frm(cs[[1]], family = cs[[2]], data = d),
                          c("did not report convergence", "boundary"))
    dp <- osa_dp(fit)
    y <- d[[cs[[3]]]]
    ok <- osa_in_step(osa_r(fit), cs[[4]](y - 1, dp), cs[[4]](y, dp))
    expect_identical(sum(!ok), 0L, info = cs[[3]])
  }
})

test_that("the residuals leave the session's random stream alone", {
  z <- osa_design()
  d <- data.frame(x = z$x, y = stats::rpois(z$n, 3))
  fit <- frm(y ~ x, family = poisson(), data = d)
  set.seed(5)
  s <- .Random.seed
  r1 <- osa_r(fit)
  expect_identical(.Random.seed, s)
  # the default uniforms are TMB's set.seed(123) ones, as the help says:
  # each residual is qnorm(F(y) - U p(y)) with U = runif(n) after
  # set.seed(123), drawn over all rows whatever the blocks (150 rows run
  # in 8 blocks here)
  expect_identical(r1, osa_r(fit, seed = 123))
  m <- osa_dp(fit)$mu
  U <- local({ set.seed(123); stats::runif(z$n) })
  ref <- stats::qnorm(stats::ppois(d$y, m) - U * stats::dpois(d$y, m))
  expect_lt(max(abs(r1 - ref)) / max(abs(ref)), 1e-10)
  # and seed = reaches them
  U7 <- local({ set.seed(7); stats::runif(z$n) })
  ref7 <- stats::qnorm(stats::ppois(d$y, m) - U7 * stats::dpois(d$y, m))
  expect_lt(max(abs(osa_r(fit, seed = 7) - ref7)) / max(abs(ref7)), 1e-10)
})

test_that("a weight between 0 and 1 gives NA and one classed warning", {
  # PIN: TMB's own qnorm(Fx - U px) warned "NaNs produced" on these rows
  # before frmtmb set them NA (dev/reviews/2026-11-06-osa.md, short
  # check, item 7)
  z <- osa_design()
  d <- data.frame(x = z$x, y = stats::rpois(z$n, exp(0.8 + 0.5 * z$x)),
                  w = rep(c(1, 0.5), length.out = z$n))
  fit <- frm(y | weights(w) ~ x, family = poisson(), data = d)
  seen <- character(0)
  r <- withCallingHandlers(osa_r(fit), warning = function(w) {
    seen <<- c(seen, if (inherits(w, "frmtmb_osa_weight")) {
      "frmtmb_osa_weight"
    } else conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, "frmtmb_osa_weight")
  expect_true(all(is.na(r[d$w < 1])))
  expect_true(all(is.finite(r[d$w == 1])))
})

test_that("a quadrature fit and the compat rows say what is true", {
  set.seed(7)
  d <- data.frame(g = factor(rep(1:10, each = 6)), x = stats::rnorm(60))
  d$y <- stats::rpois(60, exp(0.3 + 0.3 * d$x +
                                stats::rnorm(10, 0, 0.5)[d$g]))
  fq <- frm(y ~ x + (1 | g), family = poisson(), data = d,
            quadrature = TRUE)
  # PIN: 0.79.1 stopped inside TMB, "Names to be 'integrate'd must be
  # among the random parameters"
  expect_error(osa_r(fq), class = "frmtmb_osa_quadrature")
  # PIN: the row read refused while the route worked
  expect_identical(frm_compat("von_mises", "residuals_osa")$status, "works")
  expect_identical(frm_compat("tweedie", "residuals_osa")$status, "refused")
})

test_that("a case weight does not reshape a row's own step", {
  # PIN: on 0.79.1 a weight reshaped a row's own one-step distribution
  # (42 of 150 poisson rows fell outside their own step at weights of 3,
  # dev/osa-edge.R), and a weight on the "cdf" route's density term makes
  # the row's mass p^w. Where a residual sits in its step,
  # (F(y) - pnorm(r)) / p(y), is the uniform TMB drew for the row, and
  # both fits draw the same ones, whatever their estimates.
  z <- osa_design()
  d <- data.frame(x = z$x, y = stats::rpois(z$n, exp(0.8 + 0.5 * z$x)),
                  w = rep(c(1, 3), length.out = z$n))
  pos <- function(fit) {
    m <- osa_dp(fit)$mu
    (stats::ppois(d$y, m) - stats::pnorm(osa_r(fit))) /
      stats::dpois(d$y, m)
  }
  u1 <- pos(frm(y ~ x, family = poisson(), data = d))
  u2 <- pos(frm(y | weights(w) ~ x, family = poisson(), data = d))
  expect_lt(max(abs(u1 - u2)) / max(u1), 1e-6)
})

test_that("bernoulli residuals are calibrated under group-level terms", {
  # PIN: on 0.79.1 (and 0.77.0) every bernoulli residual sat near -3.7
  # on a latent fit and near -4.3 under (1 | g), seed 1
  # (dev/reviews/2026-11-02-latres.md, m6). The bound is on the mean's
  # own standard error and on the spread of a standard normal.
  set.seed(1)
  d <- data.frame(x = stats::rnorm(200), g = factor(rep(1:20, each = 10)))
  d$t <- rep(1:10, 20)
  d$bin <- stats::rbinom(200, 1, stats::plogis(0.2 + 0.4 * d$x +
                                                 stats::rnorm(20, 0, 0.6)[d$g]))
  forms <- list(bin ~ x, bin ~ x + (1 | g), bin ~ x + ar(t, g, cov = TRUE))
  for (form in forms) {
    fit <- allow_warnings(frm(form, family = bernoulli(), data = d),
                          c("did not converge", "boundary",
                            "did not report convergence"))
    r <- osa_r(fit)
    expect_true(all(is.finite(r)), info = deparse(form))
    expect_lt(abs(mean(r)) / (stats::sd(r) / sqrt(length(r))), 4)
    expect_gt(stats::sd(r), 0.75)
    expect_lt(stats::sd(r), 1.25)
  }
})

test_that("the routes that cannot read a count refuse it by name", {
  z <- osa_design()
  d <- data.frame(y = stats::rbinom(z$n, 10, z$p), x = z$x, k = 10L)
  fit <- frm(y | trials(k) ~ x, family = binomial(), data = d)
  # PIN: 0.79.1 returned NaN on 3 of 150 rows here, with no warning
  expect_error(osa_r(fit, osa_method = "oneStepGeneric"), "discreteSupport")
  # the guard's absent case: with a support the route sums the mass
  r <- osa_r(fit, osa_method = "oneStepGeneric", discreteSupport = 0:10)
  p <- osa_dp(fit)$mu
  expect_identical(sum(!osa_in_step(r, stats::pbinom(d$y - 1, 10, p),
                                    stats::pbinom(d$y, 10, p))), 0L)
  # PIN: 0.79.1 gave a gaussian approximation of a count, silently
  expect_error(osa_r(fit, osa_method = "fullGaussian"),
               "not for discrete observations")
  dt <- data.frame(y = c(rep(0, 20), stats::rgamma(60, 2, 2)), x = 1)
  ft <- frm(y ~ 1, family = tweedie(), data = dt)
  expect_error(osa_r(ft), "point mass", fixed = TRUE)
})

test_that("a non-finite residual is reported, not passed on silently", {
  z <- osa_design()
  d <- data.frame(y = stats::rpois(z$n, 3), x = z$x)
  fit <- frm(y ~ x, family = poisson(), data = d)
  # what oneStepPredict() hands back when it cannot form a row's CDF,
  # for the rows of the block it is asked for
  local_mocked_bindings(
    oneStepPredict = function(obj, ..., subset = seq_len(z$n)) {
      data.frame(residual = ifelse(subset == 1L, NaN, 0.1))
    }, .package = "RTMB")
  allow_warnings(r <- osa_r(fit), "not finite",
                 require = "1 of 150 one-step-ahead residuals are not finite")
  expect_identical(sum(is.nan(r)), 1L)
})
