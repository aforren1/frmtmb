# quantile_elf(): the extended log-F density of Fasiolo et al. (2021).
# The response is mu + lambda sigma logit(Z) with Z ~ Beta(a, b),
# a = lambda (1 - p), b = lambda p, so every reference below is formed
# from stats::dbeta(), pbeta() and qbeta() through that map, never from
# the closed form the family itself evaluates.

elf_ns <- asNamespace("frmtmb")

# the log density through the logit-beta map: log dbeta(x) plus the log
# Jacobian x (1 - x) / (lambda sigma), at x = plogis(t)
elf_ref_lpdf <- function(y, mu, sigma, p, lam) {
  a <- lam * (1 - p)
  b <- lam * p
  s <- lam * sigma
  t <- (y - mu) / s
  # dbeta() of the side below one half, since 1 - x rounds above it
  ifelse(t <= 0, stats::dbeta(stats::plogis(t), a, b, log = TRUE),
         stats::dbeta(stats::plogis(-t), b, a, log = TRUE)) +
    stats::plogis(t, log.p = TRUE) + stats::plogis(-t, log.p = TRUE) - log(s)
}

# both log tails, each from the side where its argument is below one
# half, so that plogis() never rounds to one
elf_ref_cdf <- function(y, mu, sigma, p, lam) {
  a <- lam * (1 - p)
  b <- lam * p
  t <- (y - mu) / (lam * sigma)
  lo <- stats::plogis(-abs(t))
  list(lF = ifelse(t <= 0, stats::pbeta(lo, a, b, log.p = TRUE),
                   stats::pbeta(lo, b, a, lower.tail = FALSE, log.p = TRUE)),
       lS = ifelse(t >= 0, stats::pbeta(lo, b, a, log.p = TRUE),
                   stats::pbeta(lo, a, b, lower.tail = FALSE, log.p = TRUE)))
}

elf_fam_lpdf <- function(fam, y, mu, sigma, p) {
  fam$lpdf(y, list(mu = mu, sigma = sigma, quantile = p), list())
}

test_that("quantile_elf() has asym_laplace()'s parameters and links", {
  f <- quantile_elf()
  expect_identical(f$family, "quantile_elf")
  expect_identical(f$dpars, c("mu", "sigma", "quantile"))
  expect_identical(c(f$link, f$link_sigma, f$link_quantile),
                   c("identity", "log", "logit"))
  expect_identical(quantile_elf(log)$link, "log")
  expect_error(quantile_elf("logit"),
               "not a supported link for family 'quantile_elf'")
  expect_error(quantile_elf(link_quantile = "log"),
               "not a supported link for parameter 'quantile'")
  # every spelling of the name reaches it
  expect_identical(brmsfamily("quantile_elf")$family, "quantile_elf")
  expect_identical(as_frmtmb_family("quantile_elf")$family, "quantile_elf")
})

test_that("the smoothing is one finite number above zero", {
  for (bad in list(0, -1, NA_real_, Inf, c(0.1, 0.2), "a", NULL)) {
    expect_error(quantile_elf(smoothing = bad),
                 "quantile_elf(smoothing =): the smoothing must be one",
                 fixed = TRUE, label = deparse(bad))
  }
  expect_s3_class(quantile_elf(smoothing = 1e-3), "frmtmb_family")
})

test_that("the density is the logit-beta law's and integrates to one", {
  for (lam in c(0.05, 0.25, 1, 5)) {
    for (p in c(0.1, 0.5, 0.9)) {
      fam <- quantile_elf(smoothing = lam)
      y <- 0.3 + 0.7 * lam * seq(-25, 25, by = 2.5)
      expect_equal(elf_fam_lpdf(fam, y, 0.3, 0.7, p),
                   elf_ref_lpdf(y, 0.3, 0.7, p, lam), tolerance = 1e-12,
                   label = paste(lam, p))
      tot <- stats::integrate(function(z) {
        exp(elf_fam_lpdf(fam, z, 0.3, 0.7, p))
      }, -Inf, Inf, rel.tol = 1e-10)$value
      expect_equal(tot, 1, tolerance = 1e-8, label = paste(lam, p))
    }
  }
})

test_that("as the smoothing goes to zero the density is asym_laplace()'s", {
  y <- c(-3, -0.4, 0.2, 0.9, 4)
  ald <- elf_fam_lpdf(asym_laplace(), y, 0.1, 0.8, 0.3)
  gap <- vapply(c(1e-2, 1e-3, 1e-4), function(lam) {
    max(abs(elf_fam_lpdf(quantile_elf(smoothing = lam), y, 0.1, 0.8, 0.3) -
              ald))
  }, 0)
  # the normalizer differs by O(lambda) and the smooth correction by
  # lambda exp(-|u| / lambda)
  expect_lt(gap[2], gap[1] / 5)
  expect_lt(gap[3], gap[2] / 5)
})

test_that("the log CDF pair is the incomplete beta's in both tails", {
  for (lam in c(0.05, 0.25, 1, 5)) {
    for (p in c(0.05, 0.5, 0.95)) {
      fam <- quantile_elf(smoothing = lam)
      t <- c(-600, -60, -6, -0.5, 0, 0.5, 6, 60, 600)
      y <- 0.3 + lam * 0.7 * t
      got <- fam$cdf_log(y, list(mu = 0.3, sigma = 0.7, quantile = p),
                         list())
      ref <- elf_ref_cdf(y, 0.3, 0.7, p, lam)
      expect_equal(got$lF, ref$lF, tolerance = 1e-10, label = paste(lam, p))
      expect_equal(got$lS, ref$lS, tolerance = 1e-10, label = paste(lam, p))
    }
  }
})

test_that("past 690 kernel widths each tail falls by its exact slope", {
  # I_x(a, b) = x^a / (a B(a, b)) (1 + O(x)) and log x = t + O(e^t), so
  # far below mu log F is a t - log(a B(a, b)) to a double's precision
  lam <- 0.25
  p <- 0.9
  a <- lam * (1 - p)
  b <- lam * p
  dp <- list(mu = 0, sigma = 1, quantile = p)
  fam <- quantile_elf(smoothing = lam)
  t <- c(-5000, -1000, -700)
  expect_equal(fam$cdf_log(lam * t, dp, list())$lF,
               a * t - log(a) - lbeta(a, b), tolerance = 1e-13)
  t <- c(700, 1000, 5000)
  expect_equal(fam$cdf_log(lam * t, dp, list())$lS,
               -b * t - log(b) - lbeta(a, b), tolerance = 1e-13)
  # and the two continuations meet the pair read inside without a step
  t <- 690 + c(-1e-6, 1e-6)
  ls <- fam$cdf_log(lam * t, dp, list())$lS
  expect_equal(diff(ls), -b * 2e-6, tolerance = 1e-6)
})

test_that("the taped log CDF pair has the derivatives -f / F and f / S", {
  # PIN: the continuation past 690 widths once reached the derivative
  # inside it, where F is small: d lF / d mu was off by about eps / F,
  # and exactly 0 or doubled below F = 1e-16 (seen failing on the first
  # build of this family, dev/reviews/2026-11-17-squant.md, item 1)
  t <- c(-700, -690, -689, -500, -300, -200, -150, -100, -40, -10, -1, 0,
         1, 10, 40, 100, 150, 200, 300, 500, 689, 690, 700)
  worst <- 0
  for (p in c(0.1, 0.5, 0.9)) {
    for (lam in c(0.25, 1)) {
      fam <- quantile_elf(smoothing = lam)
      q <- lam * t
      dp <- function(m) list(mu = m + 0 * q, sigma = 1, quantile = p)
      for (side in c("lF", "lS")) {
        tp <- RTMB::MakeTape(function(m) fam$cdf_log(q, dp(m), list())[[side]],
                             0)
        g <- as.numeric(tp$jacobian(0))
        v <- fam$cdf_log(q, dp(0), list())[[side]]
        # d log F / d mu = -f / F and d log S / d mu = f / S, from the
        # density and the value on doubles
        an <- (if (side == "lF") -1 else 1) *
          exp(fam$lpdf(q, dp(0), list()) - v)
        expect_true(all(is.finite(g)), label = paste(p, lam, side))
        # against the steepest slope a log tail has, max(a, b) / lambda;
        # where the derivative is below eps of it (a tail within 1e-31
        # of one) the tape gives 0, and a relative error means nothing
        slope <- max(p, 1 - p)
        expect_lt(max(abs(g - an)) / slope, 1e-12,
                  label = paste(p, lam, side))
        big <- abs(an) > .Machine$double.eps * slope
        worst <- max(worst, abs(g[big] / an[big] - 1))
      }
    }
  }
  expect_lt(worst, 1e-8)
})

test_that("the quantile function inverts the CDF", {
  for (lam in c(0.05, 0.25, 2)) {
    for (p in c(0.1, 0.5, 0.9)) {
      pr <- c(1e-12, 1e-4, 0.1, 0.5, 0.9, 1 - 1e-6)
      q <- elf_ns$elf_quantile(pr, 0.3, 0.7, p, lam)
      lf <- quantile_elf(smoothing = lam)$cdf_log(
        q, list(mu = 0.3, sigma = 0.7, quantile = p), list())$lF
      expect_equal(exp(lf), pr, tolerance = 1e-8, label = paste(lam, p))
    }
  }
  # PIN: qbeta() alone gave 177.44568, 24.556743 and 0.70978271 for these
  # (dev/reviews/2026-11-17-squant.md, item 5); the roots are the
  # review's uniroot() on the log CDF, and p 0.99 is the mirror image
  roots <- c(229.34003, 68.313613, 68.309686)
  got <- c(elf_ns$elf_quantile(0.9, 0, 1, 0.01, 0.25),
           elf_ns$elf_quantile(0.5, 0, 1, 0.01, 0.05),
           elf_ns$elf_quantile(0.5, 0, 1, 0.01, 0.001))
  expect_equal(got, roots, tolerance = 1e-7)
  mirror <- c(elf_ns$elf_quantile(0.1, 0, 1, 0.99, 0.25),
              elf_ns$elf_quantile(0.5, 0, 1, 0.99, 0.05))
  expect_equal(mirror, -roots[1:2], tolerance = 1e-7)
})

test_that("the tape has the doubles' value and finite curvature at the kink", {
  # at y = mu the check loss has a kink; the ELF's log density does not,
  # and its second derivative in mu is -1 / (4 lambda sigma^2) there
  lam <- 0.25
  sig <- 0.7
  p <- 0.8
  fam <- quantile_elf(smoothing = lam)
  f <- function(m) {
    sum(fam$lpdf(c(1, 1.3), list(mu = m + c(0, 0), sigma = sig,
                                 quantile = p), list()))
  }
  tp <- RTMB::MakeTape(f, 1)
  expect_equal(tp(1), f(1), tolerance = 1e-14)
  g <- tp$jacfun()
  h <- g$jacfun()
  # y = 1 sits on the kink: its score is (p - 1/2) / sigma
  t2 <- (1.3 - 1) / (lam * sig)
  expect_equal(as.numeric(g(1)),
               ((p - 0.5) + (stats::plogis(t2) - (1 - p))) / sig,
               tolerance = 1e-12)
  expect_equal(as.numeric(h(1)), -(0.25 + stats::dlogis(t2)) / (lam * sig^2),
               tolerance = 1e-12)
})

test_that("a fit reaches the quantile regression as the smoothing shrinks", {
  skip_if_not_installed("quantreg")
  set.seed(41)
  n <- 300
  d <- data.frame(x = stats::runif(n))
  d$y <- 1 + 2 * d$x + (0.5 + d$x) * stats::rnorm(n)
  r <- quantreg::rq(y ~ x, tau = 0.75, data = d)
  se <- summary(r, se = "nid")$coefficients[, 2]
  gap <- vapply(c(1, 0.25, 0.02), function(lam) {
    expect_no_warning(fit <- frm(bf(y ~ x, quantile = 0.75),
                                 family = quantile_elf(smoothing = lam),
                                 data = d))
    max(abs(fixef(fit)[, "Estimate"] - stats::coef(r)) / se)
  }, 0)
  expect_lt(gap[3], gap[1])
  expect_lt(gap[3], 0.25)
})

test_that("the fit solves the smoothed quantile's estimating equation", {
  # at the maximum the mean of plogis((y - mu) / (lambda sigma)) is
  # 1 - p: mu is the p quantile of y smoothed by a logistic kernel
  set.seed(42)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 0.5 + d$x + stats::rexp(n)
  for (p in c(0.2, 0.9)) {
    expect_no_warning(fit <- frm(bf(y ~ x, quantile = p),
                                 family = quantile_elf(), data = d))
    dp <- lapply(elf_ns$eval_dpars(fit)[[1L]], as.numeric)
    t <- (d$y - dp$mu) / (0.25 * dp$sigma)
    expect_equal(mean(stats::plogis(t)), 1 - p, tolerance = 1e-4)
    # the slope's score, too
    expect_lt(abs(mean(d$x * (stats::plogis(t) - (1 - p)))),
              1e-4 * mean(abs(d$x)))
  }
})

test_that("a fit has finite Hessian and sandwich standard errors", {
  set.seed(43)
  n <- 300
  d <- data.frame(x = stats::runif(n), row = seq_len(n))
  d$y <- 1 + 2 * d$x + stats::rnorm(n)
  expect_no_warning(fit <- frm(bf(y ~ x, quantile = 0.9),
                               family = quantile_elf(), data = d))
  expect_identical(fit$opt$convergence, 0L)
  se <- sqrt(diag(vcov(fit)))
  expect_true(all(is.finite(se) & se > 0))
  ss <- sqrt(diag(vcov_cluster(fit, ~ row, "CR0")))
  expect_true(all(is.finite(ss) & ss > 0))
})

test_that("cens() rows read the incomplete beta's tails", {
  set.seed(44)
  n <- 200
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + d$x + stats::rnorm(n)
  d$c <- ifelse(d$y > 2, 1, ifelse(d$y < -0.5, -1, 0))
  d$yc <- pmin(pmax(d$y, -0.5), 2)
  expect_no_warning(fit <- frm(bf(yc | cens(c) ~ x, quantile = 0.6),
                               family = quantile_elf(), data = d))
  dp <- lapply(elf_ns$eval_dpars(fit)[[1L]], as.numeric)
  ref <- elf_ref_cdf(d$yc, dp$mu, dp$sigma, 0.6, 0.25)
  ll <- ifelse(d$c == 0, elf_ref_lpdf(d$yc, dp$mu, dp$sigma, 0.6, 0.25),
               ifelse(d$c == 1, ref$lS, ref$lF))
  expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-10)
})

test_that("fitted() is the distribution's mean and osa residuals its CDF", {
  set.seed(45)
  n <- 150
  d <- data.frame(x = stats::rnorm(n))
  d$y <- d$x + stats::rnorm(n)
  expect_no_warning(fit <- frm(bf(y ~ x, quantile = 0.7),
                               family = quantile_elf(smoothing = 0.5),
                               data = d))
  dp <- lapply(elf_ns$eval_dpars(fit)[[1L]], as.numeric)
  m <- vapply(1:3, function(i) {
    stats::integrate(function(z) {
      v <- z * exp(elf_ref_lpdf(z, dp$mu[i], dp$sigma[i], 0.7, 0.5))
      # exp() of the far tails is 0, and z there is infinite
      ifelse(is.finite(v), v, 0)
    }, -Inf, Inf, rel.tol = 1e-11)$value
  }, 0)
  expect_equal(unname(fitted(fit)[1:3, "Estimate"]), m, tolerance = 1e-7)
  r <- residuals(fit, type = "osa")[, "Estimate"]
  ref <- stats::qnorm(elf_ref_cdf(d$y, dp$mu, dp$sigma, 0.7, 0.5)$lF,
                      log.p = TRUE)
  expect_equal(unname(r), ref, tolerance = 1e-8)
  s <- simulate(fit, nsim = 2, seed = 1)
  expect_identical(dim(s), c(as.integer(n), 2L))
  expect_true(all(is.finite(as.matrix(s))))
})

test_that("with (1 | g) in mu the Laplace approximation has a maximum", {
  skip_on_cran()
  # backlog item 5449's design: 80 rows in groups of 4, where
  # asym_laplace()'s Laplace objective rises without end in log sd
  set.seed(20261114)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(seq_len(n / 4),
                                                      each = 4)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n / 4, 0, 0.7)[as.integer(d$g)] +
    stats::rnorm(n)
  expect_no_warning(fit <- frm(bf(y ~ x + (1 | g), quantile = 0.5),
                               family = quantile_elf(smoothing = 0.5),
                               data = d))
  expect_identical(fit$opt$convergence, 0L)
  par <- fit$opt$par
  it <- which(names(par) == "theta")[1L]
  prof <- vapply(c(-2, 2), function(k) {
    o <- stats::nlminb(par[-it], function(p) {
      full <- par
      full[-it] <- p
      full[it] <- par[it] + k
      fit$obj$fn(full)
    })
    -o$objective
  }, 0)
  ll <- -fit$opt$objective
  expect_true(all(prof < ll))
})

test_that("weights above one with a row-level effect are refused", {
  # the kernel width is smoothing x sigma, so sigma -> 0 concentrates
  # every row's density on its own effect: the weighted likelihood rises
  # without end at weights 1.2 and has a maximum at 0.9
  set.seed(1)
  n <- 60
  d <- data.frame(x = stats::rnorm(n), obs = factor(seq_len(n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  d$w <- 1.2
  expect_error(frm(bf(y | weights(w) ~ x + (1 | obs), quantile = 0.5),
                   family = quantile_elf(), data = d),
               "The weighted likelihood of this model has no maximum",
               fixed = TRUE)
  d$w <- 0.9
  expect_no_warning(fit <- frm(bf(y | weights(w) ~ x + (1 | obs),
                                  quantile = 0.5),
                               family = quantile_elf(), data = d))
  expect_true(is.finite(as.numeric(logLik(fit))))
})
