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

test_that("a nonlinear body keeps the ill-conditioned verdict", {
  # each coefficient is its own nonlinear parameter, so there is no
  # design to take coordinates from; the SE check's "illcond" reason
  # and the NaN stay (dev/illse-nl.R)
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  f <- bf(y ~ b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 + b5 * x^5 +
            b6 * x^6 + b7 * x^7 + b8 * x^8,
          b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 ~ 1, nl = TRUE)
  r <- reparam_conds(frm(f, data = d))
  expect_true(any(grepl("the design is ill-conditioned", r$warnings,
                        fixed = TRUE)))
  lost <- sdr_of(r$value)$se_lost
  expect_length(lost, 9L)
  expect_true(all(lost == "illcond"))
  expect_null(sdr_of(r$value)$cov_fixed_root)
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
