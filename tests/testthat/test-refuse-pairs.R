# Pairs brms 2.23.0 builds that frm() refused with no recorded reason, and
# the ones it still refuses, each with its reason. The L1 identities are
# in dev/refuse-findings.md; these hold what a fit reads against an
# independent reference, and construct the case where the guard is absent.

rp_data <- function(seed = 4040, n = 160) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                  g = factor(rep(sprintf("g%02d", 1:8), length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(n, 0, 0.6)
  d$s <- 0.2 + 0.3 * abs(d$z)
  d
}

test_that("se() on skew_normal is the known sd in the density", {
  dd <- rp_data()
  fit <- frm(bf(y | se(s) ~ x), family = skew_normal(), data = dd)
  fe <- frmtmb:::fixef_by_dpar(fit)
  # sigma is mapped out, as under gaussian: mu's two and alpha are left
  expect_length(fit$opt$par, 3L)
  a <- fe$alpha[[1L]]
  mu <- fe$mu[[1L]] + fe$mu[[2L]] * dd$x
  ref <- sum(RTMBdist::dskewnorm2(dd$y, mu, dd$s, a, log = TRUE))
  expect_lt(abs(as.numeric(logLik(fit)) - ref) / abs(ref), 1e-8)
  # sigma = TRUE adds the two in quadrature, and sigma stays a parameter
  f2 <- frm(bf(y | se(s, sigma = TRUE) ~ x), family = skew_normal(),
            data = dd)
  fe2 <- frmtmb:::fixef_by_dpar(f2)
  sg <- exp(fe2$sigma[[1L]])
  mu2 <- fe2$mu[[1L]] + fe2$mu[[2L]] * dd$x
  ref2 <- sum(RTMBdist::dskewnorm2(dd$y, mu2, sqrt(sg^2 + dd$s^2),
                                   fe2$alpha[[1L]], log = TRUE))
  expect_lt(abs(as.numeric(logLik(f2)) - ref2) / abs(ref2), 1e-8)
  # another family that does not declare it still says so
  expect_error(frm(bf(y | se(s) ~ x), family = exgaussian(), data = dd,
                   dry_run = "frame"), "se()", fixed = TRUE)
})

test_that("mi() with se() reads the latent value under the known sd", {
  dd <- rp_data()
  dm <- dd
  dm$y[seq(5, 155, by = 15)] <- NA
  fit <- frm(bf(y | mi() + se(s) ~ x), data = dm)
  # a gaussian's Laplace integral is exact, so the missing rows are absent
  ref <- frm(bf(y | se(s) ~ x), data = dd[!is.na(dm$y), ])
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))) /
              abs(as.numeric(logLik(ref))), 1e-8)
  # the other two still refuse, for their reasons
  dm$ce <- 0L
  expect_error(frm(bf(y | mi() + cens(ce) ~ x), data = dm, dry_run = "frame"),
               "placeholder", fixed = TRUE)
  expect_error(frm(bf(y | mi() + trunc(lb = -50) ~ x), data = dm,
                   dry_run = "frame"), "truncation window", fixed = TRUE)
})

test_that("gr(id =, cor = FALSE) is the term with no key", {
  dd <- rp_data()
  a <- frm(bf(y ~ x + (1 + w | gr(g, id = "p", cor = FALSE)),
              sigma ~ (1 | gr(g, id = "p", cor = FALSE))), data = dd)
  b <- frm(bf(y ~ x + (1 + w || g), sigma ~ (1 | g)), data = dd)
  expect_identical(as.numeric(logLik(a)), as.numeric(logLik(b)))
  # a key with correlation still merges, so the drop above is not general
  m <- frm(bf(y ~ x + (1 | p | g), sigma ~ (1 | p | g)), data = dd,
           dry_run = "frame")
  expect_length(m$re_blocks, 1L)
})

test_that("conflicting cor under one id is refused, in both orders", {
  # brms reads the first term's flag: the two orders are two models
  dd <- rp_data()
  f <- function(form) frm(form, data = dd, dry_run = "frame")
  expect_error(f(bf(y ~ x + (1 | gr(g, id = "p")),
                    sigma ~ (1 | gr(g, id = "p", cor = FALSE)))),
               "will not pick one silently", fixed = TRUE)
  expect_error(f(bf(y ~ x + (1 | gr(g, id = "p", cor = FALSE)),
                    sigma ~ (1 | gr(g, id = "p")))),
               "will not pick one silently", fixed = TRUE)
  expect_error(f(bf(y ~ x + (1 | p | gr(g, cor = FALSE)),
                    sigma ~ (1 | p | g))),
               "will not pick one silently", fixed = TRUE)
  # another key, or another factor, does not conflict
  dd$h <- factor(rep(1:4, length.out = nrow(dd)))
  fr <- frm(bf(y ~ x + (1 | gr(g, id = "p")),
               sigma ~ (1 | gr(g, id = "q", cor = FALSE))),
            data = dd, dry_run = "frame")
  expect_length(fr$re_blocks, 2L)
})

test_that("cens() on a mixture is the mixture of the tails", {
  set.seed(7)
  n <- 200
  d <- data.frame(x = rnorm(n))
  k <- rbinom(n, 1, 0.4)
  d$y <- ifelse(k == 1, 3 + 0.5 * d$x, -3 - 0.5 * d$x) + rnorm(n)
  d$ce <- 0L
  d$ce[d$y > 3.5] <- 1L
  d$ce[d$y < -3.5] <- -1L
  fit <- frm(bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x),
             family = mixture(gaussian(), gaussian()), data = d)
  fe <- frmtmb:::fixef_by_dpar(fit)
  m1 <- fe$mu1[[1L]] + fe$mu1[[2L]] * d$x
  m2 <- fe$mu2[[1L]] + fe$mu2[[2L]] * d$x
  s1 <- exp(fe$sigma1[[1L]])
  s2 <- exp(fe$sigma2[[1L]])
  p1 <- stats::plogis(fe$theta1[[1L]])
  lse <- function(a, b) pmax(a, b) + log1p(exp(pmin(a, b) - pmax(a, b)))
  ld <- function(f, ...) lse(log(p1) + f(d$y, m1, s1, ...),
                             log(1 - p1) + f(d$y, m2, s2, ...))
  ref <- ifelse(d$ce == 0L, ld(stats::dnorm, log = TRUE),
                ifelse(d$ce == 1L,
                       ld(stats::pnorm, lower.tail = FALSE, log.p = TRUE),
                       ld(stats::pnorm, log.p = TRUE)))
  expect_gt(sum(d$ce != 0L), 0L)
  expect_lt(abs(as.numeric(logLik(fit)) - sum(ref)) / abs(sum(ref)), 1e-8)
  # truncation refuses for its reason, and a group-level mixture refuses
  # censoring
  expect_error(frm(bf(y | trunc(lb = -50) ~ 1, mu1 ~ x, mu2 ~ x),
                   family = mixture(gaussian(), gaussian()), data = d,
                   dry_run = "frame"), "truncates each component",
               fixed = TRUE)
  d$g <- factor(rep(1:20, 10))
  expect_error(frm(bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x),
                   family = mixture(gaussian(), gaussian(), groups = ~g),
                   data = d, dry_run = "frame"),
               "group-level mixture", fixed = TRUE)
})

test_that("mi() on a mixture: unread rows fit, read rows are refused", {
  set.seed(7)
  n <- 200
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  k <- rbinom(n, 1, 0.4)
  dd$y <- ifelse(k == 1, 3 + 0.5 * dd$x, -3 - 0.5 * dd$x) + rnorm(n)
  dd$y2 <- dd$y
  dd$y[1:5] <- NA
  mix <- mixture(gaussian(), gaussian())
  a <- frm(bf(y | mi() ~ 1, mu1 ~ x, mu2 ~ x), family = mix, data = dd)
  b <- frm(bf(y ~ 1, mu1 ~ x, mu2 ~ x), family = mix,
           data = dd[!is.na(dd$y), ])
  expect_lt(abs(as.numeric(logLik(a)) - as.numeric(logLik(b))) /
              abs(as.numeric(logLik(b))), 1e-6)
  f <- mvbf(bf(z ~ mi(y)), bf(y | mi() ~ 1, mu1 ~ x, mu2 ~ x, family = mix),
            rescor = FALSE)
  expect_error(frm(f, data = dd), "one mode per component", fixed = TRUE)
  expect_no_error(frm(f, data = dd, dry_run = "objective"))
})

test_that("a mixture builds under mi() on a bounded support", {
  # the placeholder at a missing row is in the support, so a component's
  # valid_y reads it as a value
  set.seed(8)
  n <- 120
  dd <- data.frame(x = rnorm(n))
  dd$y <- exp(ifelse(rbinom(n, 1, 0.5) == 1, 1, -1) + rnorm(n, 0, 0.3))
  dd$p <- stats::plogis(ifelse(rbinom(n, 1, 0.5) == 1, 1.5, -1.5) +
                          rnorm(n, 0, 0.3))
  dd$y[1:5] <- NA
  dd$p[1:5] <- NA
  for (mix in list(mixture(lognormal(), lognormal()),
                   mixture(stats::Gamma(link = "log"),
                           stats::Gamma(link = "log")))) {
    fr <- frm(bf(y | mi() ~ 1), family = mix, data = dd, dry_run = "frame")
    expect_identical(fr$mi_map$y$chart, "log")
    expect_true(all(fr$y$y > 0))
    expect_no_error(frm(bf(y | mi() ~ 1), family = mix, data = dd,
                        dry_run = "objective"))
  }
  fb <- frm(bf(p | mi() ~ 1), family = mixture(Beta(), Beta()), data = dd,
            dry_run = "frame")
  expect_identical(fb$mi_map$p$chart, "logit")
})

test_that("a guard on mi() x cens() or trunc() needs a latent value", {
  dd <- rp_data()
  dd$ce <- 0L
  a <- frm(bf(y | mi() + cens(ce) ~ x), data = dd)
  b <- frm(bf(y | cens(ce) ~ x), data = dd)
  expect_identical(as.numeric(logLik(a)), as.numeric(logLik(b)))
  expect_no_error(frm(bf(y | mi() + trunc(lb = -50) ~ x), data = dd,
                      dry_run = "frame"))
})

test_that("a known covariance is symmetric and positive definite", {
  set.seed(87)
  n <- 80
  lv <- letters[1:5]
  dd <- data.frame(x = rnorm(n), y = rnorm(n),
                   g1 = factor(sample(lv, n, TRUE), levels = lv),
                   g2 = factor(sample(lv, n, TRUE), levels = lv))
  A <- 0.5^abs(outer(1:5, 1:5, "-"))
  dimnames(A) <- list(lv, lv)
  B <- A
  B[1, 2] <- 0.4
  P <- A
  P[1, 2] <- P[2, 1] <- 1.5
  for (trm in c("(1 | mm(g1, g2, cov = M))", "(1 | gr(g1, cov = M))")) {
    fm <- bf(stats::as.formula(paste("y ~ x +", trm)))
    expect_error(frm(fm, data = dd, data2 = list(M = B), dry_run = "frame"),
                 "must be symmetric", fixed = TRUE)
    expect_error(frm(fm, data = dd, data2 = list(M = P), dry_run = "frame"),
                 "must be positive definite", fixed = TRUE)
    expect_no_error(frm(fm, data = dd, data2 = list(M = A),
                        dry_run = "frame"))
  }
})

test_that("thres() takes a multivariate model", {
  dd <- rp_data()
  dd$yo <- findInterval(dd$x + rnorm(nrow(dd)), c(-1, 0, 1)) + 1L
  fit <- frm(mvbf(bf(yo | thres(3) ~ x, family = cumulative()), bf(y ~ x),
                  rescor = FALSE), data = dd)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("the methods read a known mm() covariance and a chart", {
  set.seed(87)
  n <- 120
  lv <- letters[1:6]
  dd <- data.frame(x = rnorm(n),
                   g1 = factor(sample(lv, n, TRUE), levels = lv),
                   g2 = factor(sample(lv, n, TRUE), levels = lv))
  u <- stats::setNames(rnorm(6, 0, 0.8), lv)
  dd$y <- 0.5 * dd$x + 0.5 * (u[as.character(dd$g1)] +
                                u[as.character(dd$g2)]) + rnorm(n, 0, 0.5)
  A <- 0.5^abs(outer(1:6, 1:6, "-"))
  dimnames(A) <- list(lv, lv)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = dd,
             data2 = list(A = A))
  expect_identical(names(ranef(fit)), "mmg1g2")
  expect_length(simulate(fit, nsim = 1L)[[1L]], n)
  nd <- dd[1:3, ]
  nd$g1 <- factor("new")
  expect_error(fitted(fit, newdata = nd), "New levels")
  ft <- fitted(fit, newdata = nd, allow_new_levels = TRUE)
  expect_true(all(is.finite(ft[, "Estimate"])))

  set.seed(171)
  d2 <- data.frame(x = rnorm(100))
  d2$y <- rgamma(100, 3, 3 / exp(0.5 + 0.4 * d2$x))
  d2$y[1:8] <- NA
  f2 <- frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d2)
  r <- residuals(f2)
  expect_true(all(is.na(r[1:8, "Estimate"])))
  expect_true(all(is.finite(r[-(1:8), "Estimate"])))
  sim <- simulate(f2, nsim = 1L)[[1L]]
  expect_true(all(sim > 0))
})
