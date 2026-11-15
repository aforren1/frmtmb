# brms's four gp() kernels, cov = "exp_quad", "matern52", "matern32" and
# "exponential", exact and approximate. The references are brms 2.23.0's
# formulas written out here (R/formula-gp.R: cov_gp_*() and spd_gp_*()),
# with no frmtmb code in them. The same against brms's own functions and
# its compiled program: dev/gpfield-kernels.R and test-brms-lik-gp-exact.R.

gk_kern <- function(r, l, cov) {
  s <- r / l
  switch(cov,
    exp_quad = exp(-s^2 / 2),
    exponential = exp(-s),
    matern32 = (1 + sqrt(3) * s) * exp(-sqrt(3) * s),
    matern52 = (1 + sqrt(5) * s + 5 * s^2 / 3) * exp(-sqrt(5) * s))
}

# brms's spectral densities at frequencies w (M x D), iso a scalar l
gk_spd <- function(w, sdgp, l, cov) {
  D <- ncol(w)
  if (length(l) == 1L) l <- rep(l, D)
  q <- as.vector(w^2 %*% l^2)
  switch(cov,
    exp_quad = sdgp^2 * sqrt(2 * pi)^D * prod(l) * exp(-0.5 * q),
    exponential = sdgp^2 * (2^D * pi^(D / 2) * gamma((D + 1) / 2)) /
      sqrt(pi) * prod(l) * (1 + q)^(-(D + 1) / 2),
    matern32 = sdgp^2 * (2^D * pi^(D / 2) * gamma((D + 3) / 2) * 3^1.5) /
      (0.5 * sqrt(pi)) * prod(l) * (3 + q)^(-(D + 3) / 2),
    matern52 = sdgp^2 * (2^D * pi^(D / 2) * gamma((D + 5) / 2) * 5^2.5) /
      (0.75 * sqrt(pi)) * prod(l) * (5 + q)^(-(D + 5) / 2))
}

gk_data <- function(n = 40, seed = 9101) {
  set.seed(seed)
  d <- data.frame(x1 = round(stats::runif(n, 0, 10), 2),
                  x2 = round(stats::runif(n, 0, 5), 2))
  d$y <- sin(d$x1) + 0.3 * d$x2 + 0.2 * stats::rnorm(n)
  d
}

gk_kernels <- c("exp_quad", "matern52", "matern32", "exponential")

test_that("each kernel is brms's, isotropic and anisotropic", {
  d <- gk_data()
  for (cv in gk_kernels) {
    for (iso in c(TRUE, FALSE)) {
      fit <- frm(bf(y ~ gp(x1, x2, cov = cv, iso = iso)), data = d,
                 family = gaussian())
      bk <- fit$frame$re_blocks[[1L]]
      th <- c(0.3, log(1.7), if (!iso) log(0.6))
      K <- frmtmb:::gp_corr(th, bk)
      diag(K) <- diag(K) - frmtmb:::gp_nugget
      pos <- fit$frame$linpreds[[1L]]$gps[[1L]]$positions
      d1 <- abs(outer(pos[, 1], pos[, 1], "-"))
      d2 <- abs(outer(pos[, 2], pos[, 2], "-"))
      ref <- if (iso) {
        gk_kern(sqrt(d1^2 + d2^2), exp(th[2]), cv)
      } else {
        gk_kern(d1, exp(th[2]), cv) * gk_kern(d2, exp(th[3]), cv)
      }
      expect_lt(max(abs(K - ref)), 64 * .Machine$double.eps)
    }
  }
})

test_that("each kernel's spectral density is brms's", {
  d <- gk_data()
  for (cv in gk_kernels) {
    for (iso in c(TRUE, FALSE)) {
      # only the basis is read here: a coarse k = 5 leaves the length
      # scales flat on some kernels, and the fit says so
      fit <- allow_warnings(
        frm(bf(y ~ gp(x1, x2, k = 5, cov = cv, iso = iso)), data = d,
            family = gaussian()),
        "Standard errors are not available")
      bk <- fit$frame$re_blocks[[1L]]
      th <- c(0.3, log(0.4), if (!iso) log(0.7))
      v <- as.numeric(frmtmb:::hsgp_sds(th, bk$aux_omega, bk$gp_iso,
                                        bk$gp_cov))^2
      ref <- gk_spd(bk$aux_omega, exp(th[1]), exp(th[-1]), cv)
      expect_lt(max(abs(v / ref - 1)), 1e3 * .Machine$double.eps)
    }
  }
})

test_that("a fit's log-likelihood is the marginal under its kernel", {
  # the gaussian marginal written from the kernel at the fit's estimates:
  # the exact field and sigma integrate out in closed form
  d <- gk_data()
  for (cv in gk_kernels) {
    fit <- frm(bf(y ~ gp(x1, cov = cv)), data = d, family = gaussian())
    expect_identical(fit$opt$convergence, 0L)
    th <- fit$estimates$theta
    bk <- fit$frame$re_blocks[[1L]]
    V <- exp(2 * th[1]) * (gk_kern(abs(outer(d$x1, d$x1, "-")), exp(th[2]),
                                   cv) +
                             1e-12 * outer(d$x1, d$x1, "=="))
    diag(V) <- diag(V) + exp(2 * fit$estimates$betad[1])
    R <- chol(V)
    e <- backsolve(R, d$y - fit$estimates$beta[1], transpose = TRUE)
    ll <- -sum(log(diag(R))) - 0.5 * sum(e^2) - 0.5 * nrow(d) * log(2 * pi)
    expect_lt(abs(as.numeric(logLik(fit)) - ll), 1e-8 * abs(ll))
  }
})

test_that("the field's derivative in its length scales holds for each kernel", {
  # gp_field_dtheta() against a central difference of gp_field()
  d <- gk_data()
  for (cv in setdiff(gk_kernels, "exp_quad")) {
    for (iso in c(TRUE, FALSE)) {
      fit <- frm(bf(y ~ gp(x1, x2, cov = cv, iso = iso)), data = d,
                 family = gaussian())
      bk <- fit$frame$re_blocks[[1L]]
      th <- c(0.2, log(1.3), if (!iso) log(0.9))
      set.seed(3)
      z <- stats::rnorm(bk$dim)
      J <- frmtmb:::gp_field_dtheta(z, th, bk)
      for (k in seq_along(th)[-1L]) {
        h <- 1e-5
        tp <- th; tp[k] <- tp[k] + h
        tm <- th; tm[k] <- tm[k] - h
        fd <- (frmtmb:::gp_field(z, tp, bk) - frmtmb:::gp_field(z, tm, bk)) /
          (2 * h)
        # the factor's derivative times sd: J holds L u, gp_field() sd L z
        an <- exp(th[1]) * as.vector(frmtmb:::gp_factor(th, bk) %*%
                                       J[, k])
        expect_lt(max(abs(an - fd)), 1e-5 * max(abs(fd)))
      }
    }
  }
})

test_that("a Matern field predicts at its positions and between them", {
  d <- gk_data()
  fit <- frm(bf(y ~ gp(x1, cov = "matern32")), data = d, family = gaussian())
  nd <- data.frame(x1 = c(d$x1[1:3], 2.345, 11))
  p <- frm_linpred(fit, newdata = nd)
  expect_true(all(is.finite(p)))
  expect_lt(max(abs(p[1:3] - frm_linpred(fit)[1:3])),
            1e-8 * max(abs(p[1:3])))
  # far from every position the field returns to its mean, 0
  far <- frm_linpred(fit, newdata = data.frame(x1 = 1e3))
  expect_lt(abs(far - fit$estimates$beta[1]), 1e-6 * abs(far))
})

test_that("matern12 is the exponential kernel, and an unknown one is refused", {
  d <- gk_data()
  l1 <- logLik(frm(bf(y ~ gp(x1, cov = "matern12")), data = d,
                   family = gaussian()))
  l2 <- logLik(frm(bf(y ~ gp(x1, cov = "exponential")), data = d,
                   family = gaussian()))
  expect_identical(as.numeric(l1), as.numeric(l2))
  expect_error(frm(bf(y ~ gp(x1, cov = "nope")), data = d),
               "not a valid GP covariance kernel")
})

test_that("a derivative of a Matern field is refused by name", {
  d <- gk_data()
  fit <- frm(bf(y ~ gp(x1, cov = "matern52")), data = d, family = gaussian())
  expect_error(frm_extra_cov_deriv(fit, d[1:3, ], "x1"),
               "squared exponential kernel only")
})
