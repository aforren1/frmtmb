# Student-t distributed random effects, gr(g, dist = "student").
#
# The reference for a scalar t latent is exact, because each group's
# marginal likelihood is a one-dimensional integral. So most of these
# tests compare against adaptive quadrature rather than against another
# approximation. dev/tre-feasibility.md carries the full probe.

# --------------------------------------------------------------- setup

tre_data <- function(G = 30, n = 8, nu = 5, s = 1, seed = 4,
                     contaminate = NULL) {
  set.seed(seed)
  x <- rnorm(G * n)
  g <- factor(rep(seq_len(G), each = n))
  b <- if (is.finite(nu)) s * rt(G, df = nu) else s * rnorm(G)
  if (!is.null(contaminate)) {
    b <- rnorm(G)
    b[G] <- b[G] + contaminate
  }
  data.frame(y = 1 + 0.5 * x + b[g] + rnorm(G * n), x = x, g = g,
             b_true = b[g])
}

# Gauss-Hermite nodes and weights by Golub-Welsch: the eigenvalues of
# the Jacobi matrix are the nodes and the squared first components of
# its eigenvectors give the weights. Written out rather than taken from
# statmod, which is not a declared dependency. Verified against
# statmod::gauss.quad() to 2e-14 in both at K = 21, 51 and 81.
tre_gauss_hermite <- function(n) {
  k <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- sqrt(k / 2)
  J[cbind(k + 1L, k)] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(nodes = e$values[o], weights = sqrt(pi) * e$vectors[1, o]^2)
}

# Exact per-group marginal log-likelihood of the fitted model by
# adaptive Gauss-Hermite quadrature over the scalar latent. The latent
# density is not gaussian, so the rule is the general adaptive one:
# center and scale at the conditional mode and curvature, then undo the
# gaussian weight with exp(z^2).
tre_exact_loglik <- function(dat, beta, sigma, s, nu, K = 81) {
  r <- dat$y - beta[1] - beta[2] * dat$x
  nj <- as.vector(tapply(r, dat$g, length))
  Sj <- as.vector(tapply(r, dat$g, sum))
  SSj <- as.vector(tapply(r, dat$g, function(z) sum(z^2)))
  lp <- function(b) stats::dt(b / s, df = nu, log = TRUE) - log(s)
  lpg <- function(b) -(nu + 1) * b / (s^2 * nu + b^2)
  lph <- function(b) {
    d <- s^2 * nu + b^2
    -(nu + 1) * (d - 2 * b^2) / d^2
  }
  h <- function(B) {
    -(SSj - 2 * B * Sj + nj * B^2) / (2 * sigma^2) -
      nj / 2 * log(2 * pi * sigma^2) + lp(B)
  }
  newton <- function(b) {
    for (i in 1:200) {
      hh <- pmin(-nj / sigma^2 + lph(b), -1e-8)
      st <- -(-(nj * b - Sj) / sigma^2 + lpg(b)) / hh
      st <- pmax(pmin(st, 5 * (abs(b) + 1)), -5 * (abs(b) + 1))
      b <- b + st
    }
    b
  }
  m1 <- newton(Sj / nj)
  m2 <- newton(rep(0, length(nj)))
  m <- ifelse(h(m1) >= h(m2), m1, m2)
  curv <- -(-nj / sigma^2 + lph(m))
  sdj <- 1 / sqrt(pmax(curv, 1e-10))
  gh <- tre_gauss_hermite(K)
  B <- m + outer(sdj, sqrt(2) * gh$nodes)
  L <- -(SSj - 2 * B * Sj + nj * B^2) / (2 * sigma^2) -
    nj / 2 * log(2 * pi * sigma^2) + lp(B)
  L <- L + rep(log(gh$weights) + gh$nodes^2, each = length(nj)) +
    log(sqrt(2) * sdj)
  mx <- apply(L, 1, max)
  sum(mx + log(rowSums(exp(L - mx))))
}

# --------------------------------------------------- the density itself

test_that("us_t and diag_t are the multivariate t mvtnorm knows", {
  skip_if_not_installed("mvtnorm")
  reg <- covstruct_registry
  set.seed(11)
  for (d in c(1L, 2L, 3L)) {
    for (nu in c(2.5, 3, 5, 30)) {
      nlev <- 7L
      b <- rnorm(d * nlev)
      B <- t(matrix(b, nrow = d))
      th_us <- if (d == 1L) 0.2 else
        c(seq(-0.3, 0.3, length.out = d), rep(0.35, d * (d - 1) / 2))
      blk <- list(covstruct = "us_t", dim = d, n_levels = nlev,
                  dist_nu = nu, cnms = paste0("c", seq_len(d)))
      S <- unname(reg$us_t$vcov(th_us, blk))
      expect_equal(reg$us_t$nll(b, th_us, blk),
                   sum(mvtnorm::dmvt(B, sigma = S, df = nu, log = TRUE)),
                   tolerance = 1e-10)
      th_d <- seq(-0.3, 0.3, length.out = d)
      blk2 <- list(covstruct = "diag_t", dim = d, n_levels = nlev,
                   dist_nu = nu, cnms = paste0("c", seq_len(d)))
      S2 <- unname(reg$diag_t$vcov(th_d, blk2))
      expect_equal(reg$diag_t$nll(b, th_d, blk2),
                   sum(mvtnorm::dmvt(B, sigma = S2, df = nu,
                                     log = TRUE)),
                   tolerance = 1e-10)
    }
  }
})

test_that("a diag_t block is a multivariate t, not d univariate ones", {
  # brms shares one mixing variable across a level's coefficients even
  # under cor = FALSE, so the quadratic form couples them. The product
  # of independent t's is a DIFFERENT density, and this is the test
  # that would have caught implementing that one by mistake.
  reg <- covstruct_registry
  set.seed(3)
  b <- rnorm(2 * 5)
  blk <- list(covstruct = "diag_t", dim = 2L, n_levels = 5L,
              dist_nu = 4, cnms = c("a", "b"))
  th <- c(0.1, -0.2)
  indep <- sum(stats::dt(b * rep(exp(-th), 5), df = 4, log = TRUE)) -
    5 * sum(th)
  expect_false(isTRUE(all.equal(reg$diag_t$nll(b, th, blk), indep,
                                tolerance = 1e-6)))
  # ... and at d = 1 the two agree exactly, which is why d = 1 alone
  # could not have told them apart
  blk1 <- list(covstruct = "diag_t", dim = 1L, n_levels = 10L,
               dist_nu = 4, cnms = "a")
  expect_equal(reg$diag_t$nll(b, 0.1, blk1),
               sum(stats::dt(b * exp(-0.1), df = 4, log = TRUE)) -
                 10 * 0.1, tolerance = 1e-12)
})

test_that("the t densities differentiate on the tape", {
  reg <- covstruct_registry
  set.seed(5)
  fd_grad <- function(f, p, h = 1e-5) {
    vapply(seq_along(p), function(i) {
      a <- b <- p
      a[i] <- a[i] + h
      b[i] <- b[i] - h
      (f(a) - f(b)) / (2 * h)
    }, 0)
  }
  for (cs in c("us_t", "diag_t")) {
    for (d in c(1L, 2L, 3L)) {
      np <- reg[[cs]]$npar(d)
      nlev <- 6L
      blk <- list(covstruct = cs, dim = d, n_levels = nlev,
                  dist_nu = 3, cnms = paste0("c", seq_len(d)))
      p0 <- c(rnorm(d * nlev),
              c(seq(-0.2, 0.2, length.out = d),
                rep(0.3, max(np - d, 0)))[seq_len(np)])
      f <- function(p) {
        reg[[cs]]$nll(p[seq_len(d * nlev)],
                      p[d * nlev + seq_len(np)], blk)
      }
      tp <- RTMB::MakeTape(f, p0)
      expect_lt(max(abs(tp$jacobian(p0) - fd_grad(f, p0))), 1e-6)
    }
  }
})

# ------------------------------------------------------- the grammar

test_that("gr(dist = ) picks the covariance structure and nu", {
  d <- tre_data(G = 12, n = 5)
  f <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
           family = gaussian(), data = d)
  bk <- f$frame$re_blocks[[1L]]
  expect_identical(bk$covstruct, "us_t")
  expect_identical(bk$dist_nu, student_nu_default)

  f3 <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))),
            family = gaussian(), data = d)
  expect_identical(f3$frame$re_blocks[[1L]]$dist_nu, 3)

  fd <- frm(bf(y ~ x + diag(x | gr(g, dist = "student"))),
            family = gaussian(), data = d)
  expect_identical(fd$frame$re_blocks[[1L]]$covstruct, "diag_t")

  fc <- frm(bf(y ~ x + (x | gr(g, dist = "student"))),
            family = gaussian(), data = d)
  expect_identical(fc$frame$re_blocks[[1L]]$covstruct, "us_t")
  expect_length(fc$estimates$theta, 3L)

  # brms's explicit default is a no-op, not an error
  fg <- frm(bf(y ~ x + (1 | gr(g, dist = "gaussian"))),
            family = gaussian(), data = d)
  expect_identical(fg$frame$re_blocks[[1L]]$covstruct, "us")
  expect_null(fg$frame$re_blocks[[1L]]$dist_nu)
})

test_that("gr(dist = ) refuses what has no closed-form density", {
  d <- tre_data(G = 12, n = 5)
  A <- diag(12)
  dimnames(A) <- list(levels(d$g), levels(d$g))
  expect_error(frm(bf(y ~ (1 | gr(g, dist = "poisson"))),
                   family = gaussian(), data = d),
               "gaussian.*student")
  expect_error(frm(bf(y ~ (1 | gr(g, dist = "student", dist_nu = 2))),
                   family = gaussian(), data = d),
               "must exceed 2")
  expect_error(frm(bf(y ~ (1 | gr(g, dist_nu = 4))),
                   family = gaussian(), data = d),
               "only means something next to")
  expect_error(frm(bf(y ~ (1 | gr(g, dist = "gaussian", dist_nu = 4))),
                   family = gaussian(), data = d),
               "which has none")
  expect_error(frm(bf(y ~ (1 | gr(g, cov = A, dist = "student"))),
                   family = gaussian(), data = d, data2 = list(A = A)),
               "not a multivariate t")
  expect_error(frm(bf(y ~ cs(x | gr(g, dist = "student"))),
                   family = gaussian(), data = d),
               "us. and diag structures only")
  expect_error(frm(bf(y ~ (1 | q | gr(g, dist = "student")) +
                        (0 + x | q | gr(g, dist = "student"))),
                   family = gaussian(), data = d),
               "cannot take dist")
  d2 <- d
  d2$g2 <- factor(rep(1:6, 10))
  expect_error(frm(bf(y ~ (1 | gr(mm(g, g2), dist = "student"))),
                   family = gaussian(), data = d2),
               "loads several levels at once")
  # us_t is not a formula spelling: it is reached only through gr()
  expect_false("us_t" %in% names(formals(frm)))
})

# ----------------------------------------------- against the reference

test_that("the fitted objective is the Laplace one, to the reference", {
  d <- tre_data(G = 30, n = 8, nu = 5)
  f <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
           family = gaussian(), data = d)
  beta <- unname(f$estimates$beta)
  sigma <- exp(unname(f$estimates$betad))
  s <- exp(unname(f$estimates$theta))
  ex <- tre_exact_loglik(d, beta, sigma, s, 5)
  # 30 groups of 8: probe A1 puts the Laplace error at about -0.02 for
  # this shape, and the sign is always the same
  expect_lt(as.numeric(logLik(f)) - ex, 0)
  expect_lt(abs(as.numeric(logLik(f)) - ex), 0.1)
})

test_that("quadrature = TRUE marginalizes a t latent exactly", {
  d <- tre_data(G = 30, n = 5, nu = 3)
  fq <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))),
            family = gaussian(), data = d, quadrature = TRUE)
  beta <- unname(fq$estimates$beta)
  sigma <- exp(unname(fq$estimates$betad))
  s <- exp(unname(fq$estimates$theta))
  expect_equal(as.numeric(logLik(fq)),
               tre_exact_loglik(d, beta, sigma, s, 3),
               tolerance = 1e-6)
  # and the Laplace fit sits below it, by the amount the memo predicts
  fl <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))),
            family = gaussian(), data = d)
  expect_lt(exp(unname(fl$estimates$theta)) -
              exp(unname(fq$estimates$theta)), 0.2)
  expect_gt(exp(unname(fl$estimates$theta)) -
              exp(unname(fq$estimates$theta)), -1e-3)
})

test_that("a large nu reproduces the gaussian block", {
  d <- tre_data(G = 25, n = 8, nu = Inf)
  fn <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  ft <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 1e8))),
            family = gaussian(), data = d)
  expect_lt(max(abs(c(fn$estimates$beta, fn$estimates$betad,
                      fn$estimates$theta) -
                      c(ft$estimates$beta, ft$estimates$betad,
                        ft$estimates$theta))), 1e-6)
  expect_lt(abs(as.numeric(logLik(fn)) - as.numeric(logLik(ft))), 1e-5)

  # For a correlated block the estimate comparison is dominated by a
  # flat correlation direction, so the DENSITY is what is compared: it
  # is the exact statement anyway.
  reg <- covstruct_registry
  set.seed(9)
  for (dd in c(2L, 3L)) {
    nlev <- 20L
    b <- rnorm(dd * nlev)
    th <- c(seq(-0.3, 0.3, length.out = dd),
            rep(0.35, dd * (dd - 1) / 2))
    mk <- function(cs, nu) {
      list(covstruct = cs, dim = dd, n_levels = nlev, dist_nu = nu,
           cnms = paste0("c", seq_len(dd)))
    }
    expect_lt(abs(reg$us_t$nll(b, th, mk("us_t", 1e8)) -
                    reg$us$nll(b, th, mk("us", NULL))), 1e-5)
    th2 <- seq(-0.3, 0.3, length.out = dd)
    expect_lt(abs(reg$diag_t$nll(b, th2, mk("diag_t", 1e8)) -
                    reg$diag$nll(b, th2, mk("diag", NULL))), 1e-5)
  }
})

# ------------------------------------------------------ what it reports

test_that("VarCorr reports the SCALE and says so", {
  d <- tre_data(G = 20, n = 8, nu = 5)
  f <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 6))),
           family = gaussian(), data = d)
  V <- VarCorr(f)[[1L]]
  expect_identical(attr(V, "dist_nu"), 6)
  expect_equal(sqrt(V[1, 1]), exp(unname(f$estimates$theta)),
               tolerance = 1e-10)
  out <- capture.output(print(VarCorr(f)))
  expect_true(any(grepl("Scale", out)))
  expect_true(any(grepl("Student-t latent, nu = 6", out)))
  # the documented conversion, spelled out in the printed table
  sd_row <- exp(unname(f$estimates$theta)) * sqrt(6 / 4)
  expect_true(any(grepl(format(signif(sd_row, 5)), out, fixed = TRUE)))

  # a gaussian block keeps the old two-column form exactly
  fg <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  og <- capture.output(print(VarCorr(fg)))
  expect_false(any(grepl("Scale", og)))
  expect_null(attr(VarCorr(fg)[[1L]], "dist_nu"))
})

test_that("the sd_ alias names the scale, as brms's sd_ does", {
  d <- tre_data(G = 20, n = 8, nu = 5)
  f <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
           family = gaussian(), data = d)
  expect_true("sd_g__Intercept" %in% variables(f))
  h <- hypothesis(f, "sd_g__Intercept > 0")
  expect_equal(unname(h$estimate), exp(unname(f$estimates$theta)),
               tolerance = 1e-8)
  # confint_varcorr reports the same quantity, still the scale
  cv <- confint_varcorr(f)
  expect_equal(cv$estimate[cv$type == "sd"][1L],
               exp(unname(f$estimates$theta)), tolerance = 1e-8)
})

# ------------------------------------------------------ downstream draws

test_that("simulate() draws a multivariate t, not a gaussian", {
  d <- tre_data(G = 25, n = 6, nu = 4)
  f <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 4))),
           family = gaussian(), data = d)
  s2 <- unname(VarCorr(f)[[1L]])[1, 1]
  set.seed(1)
  B <- as.vector(replicate(4000, draw_b(f)))
  # a scalar t: b^2 / scale^2 is F(1, nu)
  k <- suppressWarnings(stats::ks.test(B^2 / s2, "pf", 1, 4))
  expect_gt(k$p.value, 0.01)
  # a gaussian draw of the same scale would fail that badly
  set.seed(1)
  Bg <- stats::rnorm(length(B), 0, sqrt(s2))
  expect_lt(suppressWarnings(
    stats::ks.test(Bg^2 / s2, "pf", 1, 4))$p.value, 1e-8)
})

test_that("a correlated draw shares one mixing variable per level", {
  d <- tre_data(G = 25, n = 8, nu = 5)
  f <- frm(bf(y ~ x + (x | gr(g, dist = "student", dist_nu = 5))),
           family = gaussian(), data = f_data <- d)
  S <- unname(VarCorr(f)[[1L]])
  set.seed(2)
  B <- t(replicate(4000, matrix(draw_b(f), nrow = 2)[, 1]))
  q <- rowSums((B %*% solve(S)) * B)
  # for a multivariate t, q/d is F(d, nu); d independent univariate t's
  # would not be
  expect_gt(suppressWarnings(
    stats::ks.test(q / 2, "pf", 2, 5))$p.value, 0.01)
})

test_that("a new level carries the t's variance, not its scale", {
  d <- tre_data(G = 20, n = 8, nu = 5)
  f <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 5))),
           family = gaussian(), data = d)
  nd <- data.frame(x = 0, g = factor("new", levels = "new"))
  p <- predict(f, newdata = nd, allow_new_levels = TRUE, se.fit = TRUE)
  s2 <- unname(VarCorr(f)[[1L]])[1, 1]
  fx <- vcov(f)[1, 1]
  expect_equal(unname(p$se.fit)^2, s2 * 5 / 3 + fx, tolerance = 1e-6)
})

test_that("frm_simulate sets the scale through the natural names", {
  d <- tre_data(G = 20, n = 8, nu = 5)
  form <- bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 5)))
  sim <- frm_simulate(form, data = d, family = gaussian(), nsim = 1,
                      seed = 7,
                      newparams = list(Intercept = 0, x = 0, sigma = 1,
                                       sd_g__Intercept = 1))
  expect_length(sim[[1L]], nrow(d))
  expect_false(anyNA(sim[[1L]]))
  # the natural name writes the SCALE, so the block's own inverse map is
  # log(): the same convention VarCorr() and confint() report
  f <- frm(form, family = gaussian(), data = d)
  bk <- f$frame$re_blocks[[1L]]
  expect_equal(covstruct_registry$us_t$from_natural(2, NULL, bk),
               log(2), tolerance = 1e-12)
  # a much larger scale really does spread the simulated response
  wide <- frm_simulate(form, data = d, family = gaussian(), nsim = 1,
                       seed = 7,
                       newparams = list(Intercept = 0, x = 0, sigma = 1,
                                        sd_g__Intercept = 8))
  expect_gt(stats::sd(wide[[1L]]), 2 * stats::sd(sim[[1L]]))
})

# ------------------------------------------------------------- payoff

test_that("a t latent holds the variance component against an outlier", {
  # the robustlmm motivation as a test: 29 clean gaussian groups plus
  # one displaced by 8 SDs. Loose tolerances: this is a direction, not
  # a calibration.
  d <- tre_data(G = 30, n = 8, contaminate = 8, seed = 21)
  fn <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  ft <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 5))),
            family = gaussian(), data = d)
  sd_n <- sqrt(unname(VarCorr(fn)[[1L]])[1, 1])
  sd_t <- sqrt(unname(VarCorr(ft)[[1L]])[1, 1]) * sqrt(5 / 3)
  expect_gt(sd_n, 1.4)          # the gaussian latent has to widen
  expect_lt(sd_t, sd_n)         # the t does not, as much
  expect_lt(sd_t, 1.5)
  # the clean groups' predicted effects are less distorted
  bt <- unname(ranef(ft)[[1L]][, 1])
  bn <- unname(ranef(fn)[[1L]][, 1])
  truth <- unname(tapply(d$b_true, d$g, `[`, 1))
  clean <- seq_len(29)
  expect_lt(sqrt(mean((bt[clean] - truth[clean])^2)),
            sqrt(mean((bn[clean] - truth[clean])^2)))
})

# ------------------------------------------------------------ REML etc

test_that("REML works over a t latent and moves the scale up", {
  d <- tre_data(G = 12, n = 5, nu = 5)
  fm <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
            family = gaussian(), data = d)
  fr <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
            family = gaussian(), data = d, REML = TRUE)
  expect_identical(fr$frame$re_blocks[[1L]]$covstruct, "us_t")
  expect_gt(unname(fr$estimates$theta), unname(fm$estimates$theta))
  expect_true(all(is.finite(unname(fr$estimates$theta))))
})

test_that("the post-fit surface is NA-free over a t block", {
  d <- tre_data(G = 20, n = 6, nu = 5)
  f <- frm(bf(y ~ x + (x | gr(g, dist = "student"))),
           family = gaussian(), data = d)
  expect_false(anyNA(fitted(f)))
  expect_false(anyNA(residuals(f)))
  expect_false(anyNA(predict(f)))
  expect_false(anyNA(ranef(f)[[1L]]))
  expect_false(anyNA(confint(f)))
  expect_false(anyNA(as.data.frame(VarCorr(f))$sdcor))
  expect_equal(unname(ngrps(f)), 20L)
  expect_false(anyNA(simulate(f, nsim = 2)))
  expect_s3_class(summary(f), "summary.frmtmb_fit")
})

test_that("frm_compat knows the t blocks", {
  cp <- frm_compat("us_t", "quadrature")
  expect_identical(cp$status, "conditional")
  expect_identical(frm_compat("us_t", "gr_cov")$status, "refused")
  expect_identical(frm_compat("diag_t", "|ID|")$status, "refused")
  expect_identical(frm_compat("us_t", "REML")$status, "works")
})
