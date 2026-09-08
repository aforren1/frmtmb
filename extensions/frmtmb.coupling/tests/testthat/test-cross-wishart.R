## The density, checked against things that are true rather than against
## another implementation of the same arithmetic. There is no second
## package fitting this likelihood, so the references here are the closed
## form maximum likelihood estimate, the distribution's own first moment,
## and the analytically known coherence bias at a true coherence of zero.

cp_form <- function() {
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1, coh ~ 1, phase ~ 1)
}

## One frequency, many units, drawn straight from a named spectral matrix.
cp_draw <- function(N, n, s11 = 1.4, s22 = 0.8, coh = 0.45, phase = 0.6) {
  a <- sqrt(s11); cm <- sqrt(coh * s22); b <- sqrt((1 - coh) * s22)
  out <- vapply(seq_len(N), function(i) {
    z1 <- complex(real = stats::rnorm(n, 0, sqrt(0.5)),
                  imaginary = stats::rnorm(n, 0, sqrt(0.5)))
    z2 <- complex(real = stats::rnorm(n, 0, sqrt(0.5)),
                  imaginary = stats::rnorm(n, 0, sqrt(0.5)))
    d1 <- a * z1
    d2 <- cm * complex(modulus = 1, argument = -phase) * z1 + b * z2
    cr <- sum(d1 * Conj(d2))
    c(sum(Mod(d1)^2), sum(Mod(d2)^2), Re(cr), Im(cr))
  }, numeric(4))
  data.frame(id = factor(seq_len(N)), w11 = out[1, ], w22 = out[2, ],
             w12r = out[3, ], w12i = out[4, ], n = n)
}

## The complete pooling estimate, which the joint maximum likelihood
## estimate must equal exactly: sum(W) / (N n).
cp_pooled <- function(d) {
  t11 <- sum(d$w11); t22 <- sum(d$w22)
  tr <- sum(d$w12r); ti <- sum(d$w12i)
  N <- sum(d$n)
  c(s11 = t11 / N, s22 = t22 / N,
    coh = (tr^2 + ti^2) / (t11 * t22), phase = atan2(ti, tr))
}

cp_est <- function(fit) {
  fx <- vapply(frmtmb::fixef(fit), function(z) z[["(Intercept)"]], 0)
  c(s11 = exp(fx[["mu"]]), s22 = exp(fx[["pow2"]]),
    coh = stats::plogis(fx[["coh"]]), phase = fx[["phase"]])
}

test_that("the maximum likelihood estimate is the pooled matrix over n", {
  # There is one thing about this family that is known in closed form and
  # this is it: with an intercept on every dpar the estimate must be
  # sum(W) / (N n) exactly, on all four coordinates at once. Anything
  # wrong in the density, the links or the parameterization moves it.
  set.seed(101)
  d <- cp_draw(N = 40L, n = 8L)
  fit <- frmtmb::frm(cp_form(), family = cross_wishart(), data = d)
  expect_equal(cp_est(fit), cp_pooled(d), tolerance = 1e-6)
})

test_that("the estimate is right at several truths and several n", {
  set.seed(102)
  for (cfg in list(list(n = 2L, coh = 0.2, phase = -1.2),
                   list(n = 5L, coh = 0.75, phase = 2.5),
                   list(n = 32L, coh = 0.05, phase = 0))) {
    d <- suppressWarnings(cp_draw(N = 30L, n = cfg$n, coh = cfg$coh,
                                  phase = cfg$phase))
    fit <- suppressWarnings(
      frmtmb::frm(cp_form(), family = cross_wishart(), data = d))
    expect_equal(cp_est(fit), cp_pooled(d), tolerance = 1e-5,
                 info = paste("n =", cfg$n, "coh =", cfg$coh))
  }
})

test_that("log det S needs no determinant and the trace no complex number", {
  # The two identities the density is built on, checked against complex
  # linear algebra rather than assumed.
  set.seed(103)
  for (r in 1:50) {
    s11 <- exp(stats::rnorm(1)); s22 <- exp(stats::rnorm(1))
    ch <- stats::runif(1, 0.01, 0.98); ph <- stats::runif(1, -pi, pi)
    s12 <- sqrt(ch * s11 * s22) * complex(modulus = 1, argument = ph)
    S <- matrix(c(s11 + 0i, Conj(s12), s12, s22 + 0i), 2, 2)
    ev <- Re(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
    expect_equal(sum(log(ev)), log(s11) + log(s22) + log(1 - ch))
    expect_true(min(ev) > 0)
    w11 <- exp(stats::rnorm(1)); w22 <- exp(stats::rnorm(1))
    w12 <- complex(real = stats::rnorm(1) * 0.1,
                   imaginary = stats::rnorm(1) * 0.1)
    W <- matrix(c(w11 + 0i, Conj(w12), w12, w22 + 0i), 2, 2)
    tr_ref <- Re(sum(diag(solve(S) %*% W)))
    recross <- sqrt(ch * s11 * s22) * (cos(ph) * Re(w12) + sin(ph) * Im(w12))
    tr_mine <- (s22 * w11 + s11 * w22 - 2 * recross) / (s11 * s22 * (1 - ch))
    expect_equal(tr_mine, tr_ref)
  }
})

test_that("the log likelihood matches an independent evaluation", {
  # The density written out again, from the definition, with a complex
  # inverse and an eigenvalue log determinant. Two spellings of the same
  # statement rather than one spelling twice.
  set.seed(104)
  d <- cp_draw(N = 12L, n = 6L)
  fit <- frmtmb::frm(cp_form(), family = cross_wishart(), data = d)
  e <- cp_est(fit)
  s12 <- sqrt(e[["coh"]] * e[["s11"]] * e[["s22"]]) *
    complex(modulus = 1, argument = e[["phase"]])
  S <- matrix(c(e[["s11"]] + 0i, Conj(s12), s12, e[["s22"]] + 0i), 2, 2)
  Si <- solve(S)
  lds <- sum(log(Re(eigen(S, symmetric = TRUE, only.values = TRUE)$values)))
  ll <- 0
  for (i in seq_len(nrow(d))) {
    W <- matrix(c(d$w11[i] + 0i, complex(real = d$w12r[i], imaginary = -d$w12i[i]),
                  complex(real = d$w12r[i], imaginary = d$w12i[i]),
                  d$w22[i] + 0i), 2, 2)
    ldw <- sum(log(Re(eigen(W, symmetric = TRUE, only.values = TRUE)$values)))
    n <- d$n[i]
    ll <- ll + (n - 2) * ldw - Re(sum(diag(Si %*% W))) - n * lds -
      (log(pi) + lgamma(n) + lgamma(n - 1))
  }
  expect_equal(as.numeric(stats::logLik(fit)), ll, tolerance = 1e-6)
})

test_that("the density is finite and exact far into the coherence tail", {
  # plogis(eta) is exactly 1 in double precision from about eta = 36.74,
  # so a density that formed 1 - C by subtraction read log(0) and divided
  # by zero there, and was already wrong in the third digit at eta = 30.
  # The complement is computed from the linear predictor instead, so both
  # ends are exact. These are the values that used to be NaN.
  lp <- frmtmb.coupling:::cw_lpdf
  aterms <- list(vint1 = 16, vreal1 = 10, vreal2 = 9.9, vreal3 = 0.3)
  y <- 10.5
  # off the tape there is no .eta_coh and the plain round trip is used;
  # on the tape there is, so build the dpars both ways and compare
  ref <- function(eta) {
    # an independent evaluation: log(1 - C) and 1/(1 - C) from eta
    lcmp <- -log1p(exp(eta))
    n <- aterms$vint1; s11 <- 1.2; s22 <- 1.1
    ch <- 1 / (1 + exp(-eta))
    ldW <- log(y * aterms$vreal1 - aterms$vreal2^2 - aterms$vreal3^2)
    rec <- sqrt(ch * s11 * s22) * (cos(0.4) * aterms$vreal2 +
                                     sin(0.4) * aterms$vreal3)
    tr <- (s22 * y + s11 * aterms$vreal1 - 2 * rec) * exp(-lcmp) / (s11 * s22)
    (n - 2) * ldW - tr - n * (log(s11) + log(s22) + lcmp) -
      (log(pi) + lgamma(n) + lgamma(n - 1))
  }
  for (eta in c(0, 3, 10, 20, 30, 36.74, 40)) {
    dp <- list(mu = 1.2, pow2 = 1.1, coh = stats::plogis(eta), phase = 0.4,
               .eta_coh = eta)
    v <- lp(y, dp, aterms)
    expect_true(is.finite(v), info = paste("eta", eta))
    expect_lt(abs(v - ref(eta)) / abs(ref(eta)), 1e-12,
              label = paste("relative error at eta", eta))
  }
  # and the slope is finite there too, which is what the optimizer needs.
  # A central difference rather than numDeriv, to add no dependency for
  # one line.
  at <- function(e) lp(y, list(mu = 1.2, pow2 = 1.1, coh = stats::plogis(e),
                               phase = 0.4, .eta_coh = e), aterms)
  for (eta in c(10, 30, 40)) {
    g <- (at(eta + 1e-4) - at(eta - 1e-4)) / 2e-4
    expect_true(is.finite(g), info = paste("slope at eta", eta))
    # the density falls away as coherence approaches 1, so the slope is
    # negative and growing in magnitude
    expect_lt(g, 0, label = paste("slope at eta", eta))
  }
})

test_that("the complement never comes from a subtraction", {
  cmp <- frmtmb.coupling:::cw_complement
  for (eta in c(0, 20, 40, 100, 700)) {
    z <- cmp(list(coh = stats::plogis(eta), .eta_coh = eta))
    expect_true(is.finite(z$log), info = paste("log at eta", eta))
    expect_equal(z$log, -log1p(exp(eta)), tolerance = 1e-12,
                 info = paste("eta", eta))
  }
  # above eta = 36.74 the subtraction route is -Inf and this one is not
  expect_true(is.infinite(log1p(-stats::plogis(40))))
  expect_true(is.finite(cmp(list(coh = 1, .eta_coh = 40))$log))
  # off the tape .eta_coh is absent and the plain round trip is used
  z <- cmp(list(coh = 0.25))
  expect_equal(z$log, log(0.75))
  expect_equal(z$inv, 1 / 0.75)
})

test_that("weights scale the log likelihood", {
  set.seed(105)
  d <- cp_draw(N = 20L, n = 8L)
  f0 <- frmtmb::frm(cp_form(), family = cross_wishart(), data = d)
  d$wt <- 1
  fw <- frmtmb::frm(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) + weights(wt) ~ 1,
               pow2 ~ 1, coh ~ 1, phase ~ 1),
    family = cross_wishart(), data = d)
  expect_equal(as.numeric(stats::logLik(fw)), as.numeric(stats::logLik(f0)),
               tolerance = 1e-6)
  d$wt <- 2
  f2 <- frmtmb::frm(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) + weights(wt) ~ 1,
               pow2 ~ 1, coh ~ 1, phase ~ 1),
    family = cross_wishart(), data = d)
  expect_equal(as.numeric(stats::logLik(f2)),
               2 * as.numeric(stats::logLik(f0)), tolerance = 1e-6)
  # and the estimate is unchanged, since every row was scaled alike
  expect_equal(cp_est(f2), cp_est(f0), tolerance = 1e-5)
})

test_that("a rank-one row is refused before any tape is built", {
  set.seed(106)
  d <- cp_draw(N = 10L, n = 4L)
  d1 <- d; d1$n <- 1L
  expect_error(frmtmb::frm(cp_form(), family = cross_wishart(), data = d1),
               "at least 2 degrees of freedom")
  d2 <- d; d2$w12r[3] <- sqrt(d2$w11[3] * d2$w22[3]) * 1.001
  expect_error(frmtmb::frm(cp_form(), family = cross_wishart(), data = d2),
               "not positive definite")
  d3 <- d; d3$w11[1] <- -1
  expect_error(frmtmb::frm(cp_form(), family = cross_wishart(), data = d3),
               "positive and finite in every row")
  d4 <- d; d4$w22[2] <- 0
  expect_error(frmtmb::frm(cp_form(), family = cross_wishart(), data = d4),
               "second auto-spectrum")
})

test_that("few degrees of freedom warn rather than pass in silence", {
  set.seed(107)
  d <- cp_draw(N = 20L, n = 3L)
  expect_warning(frmtmb::frm(cp_form(), family = cross_wishart(), data = d),
                 "fewer than 4 degrees of freedom")
})

test_that("a missing vreal or vint is refused, not fitted against nothing", {
  set.seed(108)
  d <- cp_draw(N = 10L, n = 6L)
  expect_error(
    frmtmb::frm(frmtmb::bf(w11 | vint(n) ~ 1, pow2 ~ 1, coh ~ 1, phase ~ 1),
                family = cross_wishart(), data = d))
  expect_error(
    frmtmb::frm(frmtmb::bf(w11 | vreal(w22, w12r, w12i) ~ 1, pow2 ~ 1,
                           coh ~ 1, phase ~ 1),
                family = cross_wishart(), data = d))
})

test_that("cens() and trunc() are refused by name", {
  set.seed(109)
  d <- cp_draw(N = 10L, n = 6L)
  d$ev <- 0
  expect_error(
    frmtmb::frm(frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) +
                             cens(ev) ~ 1,
                           pow2 ~ 1, coh ~ 1, phase ~ 1),
                family = cross_wishart(), data = d),
    "cens")
})

test_that("the naive bias this family exists to remove is really there", {
  # The claim in the package description, checked rather than asserted.
  # At a true coherence of zero the per-unit estimate has mean exactly
  # 1/n and averaging units does not touch it, while the model estimate
  # is the pooled one and its bias falls as 1/(N n).
  skip_on_cran()
  set.seed(110)
  N <- 20L; n <- 4L; R <- 600L
  naive <- numeric(R); pooled <- numeric(R)
  for (r in seq_len(R)) {
    d <- cp_draw(N, n, coh = 1e-8, phase = 0)
    naive[r] <- mean((d$w12r^2 + d$w12i^2) / (d$w11 * d$w22))
    pooled[r] <- cp_pooled(d)[["coh"]]
  }
  expect_equal(mean(naive), 1 / n, tolerance = 0.05)
  expect_equal(mean(pooled), 1 / (N * n), tolerance = 0.25)
  expect_lt(mean(pooled), mean(naive) / 10)
})
