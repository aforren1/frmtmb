# Standard errors on an ill-conditioned fixed design: the outer Hessian
# in the coordinates of the design's QR factor, the covariance mapped
# back, and every consumer reading it through its factor
# (R/se-reparam.R). Lane illse, dev/illse-findings.md.

reparam_conds <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

# a degree-8 raw polynomial on [1, 2]: the Hessian at unit diagonal has
# a condition number near 4e16, above 1 / machine epsilon
reparam_poly_data <- function(seed, re = FALSE) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  eta <- sin(3 * d$x)
  if (re) eta <- eta + stats::rnorm(20, 0, 0.5)[d$g]
  d$y <- eta + stats::rnorm(200, 0, 0.3)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  d
}
reparam_fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
reparam_nd <- function(d) {
  nd <- data.frame(x = c(1, 1.25, 1.5, 1.75, 2), g = d$g[1])
  for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
  nd
}
reparam_products <- "a product such as a' V a"

test_that("a degree-8 raw polynomial gets lm()'s covariance and SEs", {
  # 0.70.0 lost all nine coefficients as "illcond" and every prediction
  # that read them; lm() reports them through its QR. In the design's
  # coordinates the covariance of every entry is within 3e-6 of the
  # exact one at the fit's own point (256-bit, dev/illse-check.R), and
  # the Newton steps there reach lm()'s maximum
  d <- reparam_poly_data(1)
  r <- reparam_conds(frm(reparam_fo, data = d))
  fit <- r$value
  expect_length(r$warnings, 0L)
  expect_length(sdr_of(fit)$se_lost, 0L)
  lmf <- stats::lm(reparam_fo, data = d)
  n <- nrow(d)
  k <- length(stats::coef(lmf))
  # vcov() hands out the whole matrix, so it says once what a product
  # of it loses
  V <- allow_warnings(vcov(fit), reparam_products,
                      require = reparam_products)
  Vl <- stats::vcov(lmf) * (n - k) / n
  expect_lt(max(abs(V / unname(Vl) - 1)), 1e-4)
  # the estimates are lm()'s to a small fraction of an SE
  expect_lt(max(abs(fixef(fit)[, "Estimate"] - stats::coef(lmf)) /
                  sqrt(diag(Vl))), 0.01)
  # fitted(): through the factor, against predict.lm(se.fit = TRUE)
  nd <- reparam_nd(d)
  p <- reparam_conds(fitted(fit, newdata = nd, re_formula = NA))
  expect_length(p$warnings, 0L)
  pl <- stats::predict(lmf, newdata = nd, se.fit = TRUE)
  ref <- unname(pl$se.fit) * sqrt((n - k) / n)
  se <- unname(p$value[, "Est.Error"])
  expect_lt(max(abs(se / ref - 1)), 1e-4)
  # hypothesis() of the polynomial at x = 1.5 reads the same factor
  h <- reparam_conds(hypothesis(
    fit, paste0("Intercept + ", paste0(1.5^(1:8), " * x", 1:8,
                                       collapse = " + "), " = 0"),
    class = "b"))
  expect_length(h$warnings, 0L)
  expect_lt(abs(h$value$hypothesis$Est.Error / se[3] - 1), 1e-4)
  # the parameter draws of predict() come from the factor: chol() of
  # the whole matrix has no digits left. 4000 draws put the Monte Carlo
  # error of an SD near 1.1 percent; the bound is nine of those
  set.seed(11)
  dr <- frmtmb:::predict_par_drawer(fit, TRUE, 4000L)
  delta <- vapply(seq_len(4000L), attr(dr, "delta"), numeric(10))
  A <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, nd)
  sd_draw <- apply(A %*% delta[1:9, ], 1L, stats::sd)
  expect_lt(max(abs(sd_draw / ref - 1)), 0.1)
})

test_that("emmeans and marginaleffects on an ill-conditioned design", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("marginaleffects")
  d <- reparam_poly_data(1)
  fit <- frm(reparam_fo, data = d)
  nd <- reparam_nd(d)
  se <- unname(fitted(fit, newdata = nd, re_formula = NA)[, "Est.Error"])
  # emmeans takes the grid route, G V G' formed through the factor; the
  # design route handed emmeans X and V and it lost every digit (B6 of
  # the seflags review: 0.57 to 1.30 of the true SE)
  em <- reparam_conds(summary(emmeans::emmeans(
    fit, ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
    at = as.list(nd[3, paste0("x", 1:8)]))))
  expect_length(em$warnings, 0L)
  expect_lt(abs(em$value$SE / se[3] - 1), 1e-6)
  # marginaleffects forms J V J' from get_vcov() itself, which cannot be
  # given a factor: NaN with one warning, not a wrong number
  mf <- allow_warnings(
    marginaleffects::predictions(fit, newdata = nd),
    c("marginaleffects forms standard errors",
      # insight, once per session, on a model whose data it cannot see
      "could not be extracted reliably"),
    require = "marginaleffects forms standard errors")
  expect_true(all(is.nan(mf$std.error)))
})

test_that("with (1 | g) the factor reaches the group effects", {
  # the exact prediction SE from a QR of the whitened design at the
  # fit's sigma and sd, the construction of test-se-flags.R's degree-7
  # test; 0.70.0 gave NaN here for all nine coefficients
  d <- reparam_poly_data(3, re = TRUE)
  fo <- update(reparam_fo, . ~ . + (1 | g))
  r <- reparam_conds(frm(fo, data = d))
  fit <- r$value
  expect_length(r$warnings, 0L)
  expect_length(sdr_of(fit)$se_lost, 0L)
  nd <- reparam_nd(d)
  p <- reparam_conds(fitted(fit, newdata = nd, re_formula = NA))
  expect_length(p$warnings, 0L)
  on <- frmtmb:::outer_par_names(fit)
  sg <- exp(fit$opt$par[grepl("^sigma", on)])
  tg <- exp(fit$opt$par[on == "theta_1"])
  X <- stats::model.matrix(reparam_fo, d)
  X0 <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, nd)
  Z <- stats::model.matrix(~ 0 + g, d)
  L <- chol(sg^2 * diag(nrow(X)) + tg^2 * Z %*% t(Z))
  q <- qr(backsolve(L, X, transpose = TRUE), tol = 1e-30)
  A <- backsolve(qr.R(q), t(X0[, q$pivot, drop = FALSE]), transpose = TRUE)
  ref <- sqrt(colSums(A^2))
  # conditional on sigma and the sd, which the fit's covariance is not,
  # so the bound is the one test-se-flags.R gives that construction
  expect_lt(max(abs(unname(p$value[, "Est.Error"]) / ref - 1)), 0.01)
  # the same likelihood with poly(x, 8) is well conditioned; a group's
  # own prediction carries its effect's variance through the factor too,
  # where J Vf J' formed whole gave it no digits
  ref2 <- frm(y ~ poly(x, 8) + (1 | g), data = d)
  pg <- reparam_conds(fitted(fit, newdata = nd))
  expect_length(pg$warnings, 0L)
  pr <- fitted(ref2, newdata = nd, re_formula = NA)
  pgr <- fitted(ref2, newdata = nd)
  expect_lt(max(abs(p$value[, "Est.Error"] / pr[, "Est.Error"] - 1)),
            1e-3)
  expect_lt(max(abs(pg$value[, "Est.Error"] / pgr[, "Est.Error"] - 1)),
            1e-3)
  # and so does ranef(condVar = TRUE), which read sdreport()'s own
  # inverse of the natural Hessian and gave NaN (dev/illse-routes.R)
  cs <- attr(ranef(fit, condVar = TRUE)$g, "condSD")
  csr <- attr(ranef(ref2, condVar = TRUE)$g, "condSD")
  expect_lt(max(abs(cs / csr - 1)), 1e-3)
})

test_that("an ordinal fit's category SEs on an ill-conditioned design", {
  # B7 of the seflags review: P(y = 1) SEs 15925 to 255929 times the
  # poly() fit's through the category route. The same likelihood with
  # poly(x, 8) is well conditioned and is the reference
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- factor(cut(2 * sin(3 * d$x) + stats::rlogis(200),
                    c(-Inf, -1, 0, 1, Inf), labels = FALSE),
                ordered = TRUE)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  r <- reparam_conds(frm(reparam_fo, data = d, family = cumulative()))
  expect_length(r$warnings, 0L)
  ref <- frm(y ~ poly(x, 8), data = d, family = cumulative())
  nd <- data.frame(x = c(1.2, 1.5, 1.8))
  for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
  p <- reparam_conds(fitted(r$value, newdata = nd))
  expect_length(p$warnings, 0L)
  pr <- fitted(ref, newdata = nd)
  expect_lt(max(abs(p$value[, "Est.Error", ] / pr[, "Est.Error", ] - 1)),
            1e-3)
})

test_that("a well-conditioned design keeps its covariance untouched", {
  # the guard's absent case: degree 3, condition number near 3e6 at
  # unit diagonal, far below the 1e8 above which the design's
  # coordinates are used; no factor, and vcov() says nothing
  d <- reparam_poly_data(1)
  fit <- frm(y ~ x1 + x2 + x3, data = d)
  expect_null(sdr_of(fit)$cov_fixed_root)
  v <- reparam_conds(vcov(fit))
  expect_length(v$warnings, 0L)
  expect_null(frm_lp_basis(fit, newdata = reparam_nd(d))$V_root)
})

test_that("a nonlinear body takes the coordinates of its Jacobian", {
  # each coefficient is its own nonlinear parameter, so there is no
  # design; the body's Jacobian is one, and on a body linear in its
  # parameters it is the design itself. 0.71.0 stopped 0.41 below lm()'s
  # maximum and lost all nine as "illcond" (dev/illse2-p3nl.R poly)
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  f <- bf(y ~ b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 + b5 * x^5 +
            b6 * x^6 + b7 * x^7 + b8 * x^8,
          b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 ~ 1, nl = TRUE)
  r <- reparam_conds(frm(f, data = d))
  expect_length(r$warnings, 0L)
  fit <- r$value
  expect_length(sdr_of(fit)$se_lost, 0L)
  expect_false(is.null(sdr_of(fit)$cov_fixed_root))
  lmf <- stats::lm(y ~ stats::poly(x, 8, raw = TRUE), data = d)
  ll <- as.numeric(stats::logLik(lmf))
  expect_lt(abs(as.numeric(stats::logLik(fit)) / ll - 1), 1e-8)
  n <- nrow(d)
  Vl <- stats::vcov(lmf) * (n - 9) / n
  V <- allow_warnings(vcov(fit, full = TRUE), reparam_products)[1:9, 1:9]
  expect_lt(max(abs(V / unname(Vl) - 1)), 1e-4)
})

test_that("a covariate far from zero reaches lm()'s maximum and SEs", {
  # x on [1e5, 1e5 + 1]: 0.70.0 stopped at nlminb's "false convergence"
  # 157 below lm()'s log-likelihood, at its start slope, with no SE.
  # Newton steps in the design's coordinates now run from such a stop
  # too (newton_polish()), damped while the Hessian there is indefinite
  # (dev/illse-offdbg.R)
  set.seed(2)
  d <- data.frame(x = 1e5 + stats::runif(200))
  d$y <- 1 + 2 * (d$x - 1e5) + stats::rnorm(200, 0, 0.3)
  # whether nlminb itself calls the stop converged is rounding's choice;
  # either way the warning, if any, says the stop is the maximum
  fit <- allow_warnings(frm(y ~ x, data = d),
                        "so the estimates are at the maximum")
  lmf <- stats::lm(y ~ x, data = d)
  ll <- as.numeric(stats::logLik(lmf))
  expect_lt(abs(as.numeric(stats::logLik(fit)) / ll - 1), 1e-8)
  V <- allow_warnings(vcov(fit), reparam_products)
  Vl <- stats::vcov(lmf) * (200 - 2) / 200
  expect_lt(max(abs(V / unname(Vl) - 1)), 1e-4)
})

test_that("vcov_cluster() forms the sandwich through the factor", {
  # 0.70.0 gave NaN; formed as B S'S B from the whole bread the
  # sandwich was off by 2.3 to 12 on the correlation scale and a fitted
  # value's SE through it was 0 to 1511 times the exact one (review of
  # lane illse, item 1). The reference is lm()'s QR: at the maximum the
  # bread is block diagonal between the coefficients and log sigma, so
  # the coefficients' block is the sum over clusters of u u' with
  # u = R^-1 Q_g' r_g, sigma cancelling
  d <- reparam_poly_data(1)
  fit <- frm(reparam_fo, data = d)
  r <- reparam_conds(vcov_cluster(fit, ~ g, full = TRUE))
  Vc <- r$value
  expect_false(is.null(attr(Vc, "root")))
  expect_length(grep(reparam_products, r$warnings, fixed = TRUE), 1L)
  X <- stats::model.matrix(reparam_fo, d)
  q <- qr(X, tol = 1e-30)
  expect_identical(q$pivot, seq_len(ncol(X)))
  Q <- qr.Q(q)
  Rm <- qr.R(q)
  res <- d$y - drop(X %*% fixef(fit)[, "Estimate"])
  U <- vapply(split(seq_len(nrow(d)), d$g), function(ix) {
    backsolve(Rm, drop(crossprod(Q[ix, , drop = FALSE], res[ix])))
  }, numeric(ncol(X)))
  Vref <- tcrossprod(U)
  bb <- which(frmtmb:::outer_par_map(fit)$comp == "beta")
  s <- sqrt(diag(Vref))
  expect_lt(max(abs(Vc[bb, bb] - Vref) / outer(s, s)), 1e-4)
  # hypothesis() reads the factor that came with the matrix
  a <- c(1, 1.5^(1:8))
  hs <- paste0("Intercept + ", paste0(1.5^(1:8), " * x", 1:8,
                                      collapse = " + "), " = 0")
  h <- reparam_conds(hypothesis(fit, hs, class = "b", vcov = Vc))
  expect_length(h$warnings, 0L)
  ref <- sqrt(sum(drop(crossprod(a, U))^2))
  expect_lt(abs(h$value$hypothesis$Est.Error / ref - 1), 1e-4)
  # without its factor the matrix gives no SE to a combination, and
  # says why, while a single coefficient keeps its SE
  V2 <- Vc
  attr(V2, "root") <- NULL
  h2 <- allow_warnings(
    hypothesis(fit, c(hs, "x8 = 0"), class = "b", vcov = V2),
    "the supplied matrix does not carry one",
    require = "the supplied matrix does not carry one")
  expect_true(is.nan(h2$hypothesis$Est.Error[1L]))
  expect_true(is.finite(h2$hypothesis$Est.Error[2L]))
})

test_that("insight::get_varcov() is NaN where products lose their digits", {
  # ggeffects, modelbased and parameters form a' V a from it; from the
  # whole matrix at degree 8 it was 0.779 of lm()'s SE with no warning,
  # where 0.70.0 gave NaN (review of lane illse, item 2)
  skip_if_not_installed("insight")
  d <- reparam_poly_data(1)
  fit <- frm(reparam_fo, data = d)
  V <- allow_warnings(insight::get_varcov(fit),
                      "insight::get_varcov()",
                      require = "insight::get_varcov()")
  expect_true(all(is.nan(V)))
})

test_that("marginaleffects gets NaN on every fit held as a factor", {
  # a binomial with mixed scales: condition number near 4e11, below the
  # 1e12 at which a product formed from the whole matrix loses its
  # digits, but marginaleffects' numerical Jacobian, multiplied by the
  # cancellation in J V J', gave predictions() SEs 4.8 to 15.1 times the
  # exact ones with no warning (review of lane illse, item 3, seed 201)
  skip_if_not_installed("marginaleffects")
  set.seed(201)
  n <- 300
  d <- data.frame(u = 100 + stats::runif(n), z = 1e-3 * stats::rnorm(n),
                  w = 1e4 * stats::rnorm(n))
  d$u2 <- d$u^2
  d$N <- 20L
  eta <- -0.3 + 3 * (d$u - 100.5) - 2 * (d$u - 100.5)^2 + 300 * d$z +
    5e-5 * d$w
  d$y <- stats::rbinom(n, d$N, stats::plogis(eta))
  fit <- allow_warnings(frm(y | trials(N) ~ u + u2 + z + w, data = d,
                            family = binomial()),
                        "so the estimates are at the maximum")
  root <- sdr_of(fit)$cov_fixed_root
  expect_false(is.null(root))
  expect_false(isTRUE(attr(root, "lossy")))
  nd <- data.frame(u = 100.5, z = 0, w = 0, N = 20L)
  nd$u2 <- nd$u^2
  mf <- allow_warnings(
    marginaleffects::predictions(fit, newdata = nd),
    c("marginaleffects forms standard errors",
      "could not be extracted reliably"),
    require = "marginaleffects forms standard errors")
  expect_true(all(is.nan(mf$std.error)))
})

# raw polynomials with (1 | g) on data with no group effect: the sd
# runs to zero (the review of lane illse, dev/illse-rev-edge4.R)
reparam_edge_fit <- function(deg, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(300, 1, 2),
                  g = factor(rep(1:20, length.out = 300)))
  d$y <- sin(3 * d$x) + stats::rnorm(300, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                "+ (1 | g)"))
  r <- reparam_conds(frm(fo, data = d))
  r$data <- d
  r$lm <- stats::lm(stats::as.formula(paste("y ~", paste(xn,
                                                       collapse = " + "))),
                    data = d)
  r
}

# lm()'s coefficient SEs on the ML scale: the covariance of the
# coefficients with the sd at zero
reparam_lm_se <- function(r) {
  k <- length(stats::coef(r$lm))
  n <- nrow(r$data)
  unname(sqrt(diag(stats::vcov(r$lm)) * (n - k) / n))
}

test_that("an sd near zero is judged by the check, not the route", {
  # the route took it, since its unit-diagonal eigenvalue is near 1, and
  # gave SE(log sd) 4.75e5 with no warning; 0.70.0 named it (item 4)
  r <- reparam_edge_fit(5, 4)
  lost <- sdr_of(r$value)$se_lost
  expect_identical(names(lost), "theta_1")
  expect_true(any(grepl("theta_1", r$warnings, fixed = TRUE)))
})

test_that("with the sd at its boundary the coefficients keep lm()'s SEs", {
  # tier 3 removes the sd, and the rest of a degree-8 design goes to its
  # own coordinates: in natural ones tier 3 called x3 to x6 flat and
  # kept the other coefficients' SEs at 5 to 10 percent of lm()'s
  # (item 4); 0.70.0 lost all nine as "illcond"
  r <- reparam_edge_fit(8, 1)
  lost <- sdr_of(r$value)$se_lost
  expect_identical(names(lost), "theta_1")
  expect_false(is.null(sdr_of(r$value)$cov_fixed_root))
  se <- unname(sqrt(diag(allow_warnings(vcov(r$value), reparam_products))))
  expect_lt(max(abs(se / reparam_lm_se(r) - 1)), 1e-3)
})

test_that("a fit at the maximum with an sd at zero is not unconverged", {
  # the Newton steps reach lm()'s maximum to 6.4e-9, but the sd's empty
  # row left the Hessian without a Cholesky factor, so the gradient check
  # warned "may not have converged" and every SE was NaN (item 5)
  r <- reparam_edge_fit(6, 4)
  expect_length(grep("Large maximum absolute gradient", r$warnings,
                     fixed = TRUE), 0L)
  expect_length(grep("did not report convergence", r$warnings,
                     fixed = TRUE), 0L)
  se <- unname(sqrt(diag(allow_warnings(vcov(r$value), reparam_products))))
  expect_lt(max(abs(se / reparam_lm_se(r) - 1)), 1e-3)
})

test_that("an sd stopped near zero below a small interior maximum moves", {
  # degree 7, seed 4: nlminb stopped at an sd of 1.9e-6 at code 0, but
  # the maximum is at 0.033 times sigma, 1.2e-3 higher; the headroom
  # then called the fit short and every coefficient lost its SE. The
  # optimizer restarts from the larger sd the upward probe finds
  # (sd_up_escape()). The reference is the profile of the
  # log-likelihood in tau = sd / sigma: per group, I - a J whitens the
  # rows, and the residual of the QR of the whitened design is the GLS
  # residual (dev/illse2-p2.R, 256-bit, gives the same maximum)
  r <- reparam_edge_fit(7, 4)
  expect_length(r$warnings, 0L)
  d <- r$data
  X <- stats::model.matrix(r$lm)
  gi <- split(seq_len(nrow(d)), d$g)
  prof <- function(tau) {
    yt <- d$y
    Xt <- X
    ld <- 0
    for (ix in gi) {
      ng <- length(ix)
      a <- (1 - 1 / sqrt(1 + ng * tau^2)) / ng
      yt[ix] <- d$y[ix] - a * sum(d$y[ix])
      Xt[ix, ] <- X[ix, , drop = FALSE] -
        a * rep(colSums(X[ix, , drop = FALSE]), each = ng)
      ld <- ld + log(1 + ng * tau^2)
    }
    rss <- sum(qr.resid(qr(Xt, tol = 1e-12), yt)^2)
    n <- nrow(d)
    -(n / 2) * log(2 * pi * rss / n) - n / 2 - ld / 2
  }
  top <- stats::optimize(prof, c(0, 0.3), maximum = TRUE, tol = 1e-8)
  expect_gt(top$maximum, 0.01)
  ll <- as.numeric(stats::logLik(r$value))
  expect_lt((top$objective - ll) / abs(top$objective), 1e-7)
  se <- sqrt(diag(allow_warnings(vcov(r$value), reparam_products)))
  expect_true(all(is.finite(se)))
})

test_that("an active bound on an ill-conditioned design is reached", {
  # degree 6, an upper bound one SE below lm()'s estimate of x6: every
  # full Newton step crossed the bound and was refused, and 0.71.0
  # stopped 0.12 below the constrained maximum at code 0 with no warning
  # and coefficient SEs 5.8 to 36 times the right ones. The steps now
  # hold a coordinate they would push out at its bound
  # (newton_polish()), and the others' covariance is taken in the
  # coordinates of the design without x6 (dev/illse2-p4bound.R seed 1).
  # The constrained maximum is lm() with x6 as an offset at the bound
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.2)
  for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6
  l0 <- stats::lm(fo, data = d)
  ub <- stats::coef(l0)[["x6"]] - sqrt(diag(stats::vcov(l0)))[["x6"]]
  fit <- allow_warnings(
    frm(fo, data = d, prior = set_prior("", class = "b", coef = "x6",
                                        ub = ub)),
    "x6", require = "x6")
  expect_identical(names(sdr_of(fit)$se_lost), "x6")
  d$off <- ub * d$x6
  l1 <- stats::lm(y ~ x1 + x2 + x3 + x4 + x5 + offset(off), data = d)
  ll <- as.numeric(stats::logLik(l1))
  expect_lt(abs(as.numeric(stats::logLik(fit)) / ll - 1), 1e-7)
  se <- allow_warnings(fixef(fit), "x6")[c("Intercept", paste0("x", 1:5)),
                                         "Est.Error"]
  ref <- sqrt(diag(stats::vcov(l1)) * (200 - 6) / 200)
  expect_lt(max(abs(unname(se / ref) - 1)), 1e-3)
})

test_that("raising the optimizer's caps keeps the sd escape", {
  # a caller who asks for more iterations asks for more work: with
  # iter.max = eval.max = 5000 the escape was off, and degree 7 seed 4
  # stayed 1.2e-3 short with every coefficient SE NaN (review of lane
  # illse2, finding 2). A cap below the default keeps its stop
  # (test-se-conv.R)
  r <- reparam_edge_fit(7, 4)
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  w <- character()
  f2 <- withCallingHandlers(
    frm(fo, data = r$data,
        control = frmtmb_control(optCtrl = list(iter.max = 5000,
                                                eval.max = 5000))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 0L)
  expect_true(isTRUE(f2$opt$sd_escape[["gain"]] > 1e-3))
  ll <- as.numeric(stats::logLik(r$value))
  expect_lt(abs(as.numeric(stats::logLik(f2)) / ll - 1), 1e-9)
  expect_false(frmtmb:::optctrl_stops_earlier(list(iter.max = 5000,
                                                   rel.tol = 1e-12)))
  expect_true(frmtmb:::optctrl_stops_earlier(list(iter.max = 3)))
  expect_true(frmtmb:::optctrl_stops_earlier(list(rel.tol = 1e-6)))
})

reparam_aexp_fit <- function(seed, A0) {
  set.seed(seed)
  n <- 2000
  d <- data.frame(x = stats::runif(n, 1, 2), z = stats::runif(n, 0, 10))
  for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
  d$y <- sin(3 * d$x) + A0 * exp(-0.5 * d$z) + stats::rnorm(n, 0, 0.01)
  f <- bf(y ~ b0 + b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4 + b5 * x5 +
            b6 * x6 + A * exp(-k * z),
          b0 + b1 + b2 + b3 + b4 + b5 + b6 + A + k ~ 1, nl = TRUE)
  l0 <- stats::lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data = d)
  st <- list(beta = c(unname(stats::coef(l0)), A0, 0.5))
  list(d = d, fit = frm(f, data = d, start = st))
}

test_that("an identified nonlinear body passes the asymmetry guard", {
  # a degree-6 raw polynomial beside A exp(-k z), A = 1 (t(k) near
  # 460): the route takes it, and with every column extrapolated its
  # covariance is within 3.2e-7 of the exact one (review of lane
  # illse2, findings 3 and 4). Reference: the Hessian at the maximum,
  # (J'J - S) / s2 with S the residuals times the body's second
  # derivatives in (A, k), inverted through the QR of the Jacobian J
  r <- reparam_aexp_fit(1, 1)
  fit <- r$fit
  d <- r$d
  expect_false(is.null(sdr_of(fit)$cov_fixed_root))
  b <- fit$opt$par[1:9]
  s2 <- exp(2 * fit$opt$par[[10]])
  X <- as.matrix(cbind(1, d[paste0("x", 1:6)]))
  e <- exp(-b[[9]] * d$z)
  J <- cbind(X, e, -b[[8]] * d$z * e)
  res <- d$y - drop(X %*% b[1:7]) - b[[8]] * e
  S <- matrix(0, 9, 9)
  S[8, 9] <- S[9, 8] <- sum(res * (-d$z * e))
  S[9, 9] <- sum(res * b[[8]] * d$z^2 * e)
  q <- qr(J, tol = 1e-12)
  Ri <- matrix(0, 9, 9)
  Ri[q$pivot, ] <- backsolve(qr.R(q), diag(9))
  G <- crossprod(Ri, S %*% Ri)
  Fm <- Ri %*% t(chol(solve(diag(9) - G)))
  se_ref <- sqrt(s2 * rowSums(Fm^2))
  se <- sqrt(diag(allow_warnings(vcov(fit, full = TRUE),
                                 reparam_products)))[1:9]
  expect_lt(max(abs(unname(se) / se_ref - 1)), 1e-4)
})

test_that("a nearly unidentified nonlinear body is refused by the guard", {
  # the same body with A = 0.05 (t(k) near 22), seed 2: the Jacobian's
  # coordinates are not linear enough for the differences there, and
  # the Hessian in them is asymmetric by 6.6e-3 of its largest entry
  # (at most 1.7e-5 on fits the route should take). Taken anyway, its
  # covariance was 9.1e-2 from the exact one on the scale of the
  # correlation, silently; refused, the fit keeps 0.71.0's, 5.9e-4 off
  # (dev/illse2-rev-asym.R with I2_ASYM=off)
  fit <- reparam_aexp_fit(2, 0.05)$fit
  expect_null(sdr_of(fit)$cov_fixed_root)
  p <- fit$opt$par
  M <- frmtmb:::reparam_basis(fit$frame, names(p), est = fit$estimates)
  h <- frmtmb:::reparam_hessian(fit$obj, p, M)
  expect_gt(max(h$E) / max(abs(h$H)), frmtmb:::se_reparam_asym)
})
