# REML on an ill-conditioned fixed design: beta is integrated with the
# group effects, so the objective is taped in the design's coordinates
# (beta_whitener()) and the joint covariance is held as a factor
# (sdr_unwhiten()). Lane illse2, dev/illse2-findings.md, item 3.

reml_data <- function(seed, deg) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(200, 0, 0.3)
  for (k in seq_len(deg)) d[[paste0("x", k)]] <- d$x^k
  d
}

# GLS at t = sd^2 / sigma^2: per group, I - a J whitens the rows, and the
# QR of the whitened design gives the estimate and the conditional
# covariance's factor, sigma R^-1
reml_gls <- function(X, y, g, t) {
  Xt <- X
  yt <- y
  for (ix in split(seq_along(y), g)) {
    ng <- length(ix)
    a <- (1 - 1 / sqrt(1 + ng * t)) / ng
    yt[ix] <- y[ix] - a * sum(y[ix])
    Xt[ix, ] <- X[ix, , drop = FALSE] -
      a * rep(colSums(X[ix, , drop = FALSE]), each = ng)
  }
  q <- qr(Xt, tol = 1e-12)
  stopifnot(q$rank == ncol(X))
  Ri <- matrix(0, ncol(X), ncol(X))
  Ri[q$pivot, ] <- backsolve(qr.R(q), diag(ncol(X)))
  list(beta = qr.coef(q, yt), Ri = Ri)
}

test_that("REML on a degree-7 raw polynomial gets the exact covariance", {
  # 0.71.0 stopped short with a convergence warning and coefficient SEs
  # 1.5 to 1.9 times the exact ones (degree 7, seeds 1 to 3); at degree
  # 8 every SE was NaN, and lme4 stops with "Downdated VtV is not
  # positive definite" (dev/illse2-p3reml.R, 256-bit reference)
  d <- reml_data(1, 7)
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  w <- character()
  fit <- withCallingHandlers(frm(fo, data = d, REML = TRUE),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
  expect_length(w, 0L)
  expect_false(is.null(sdr_of(fit)$joint_root))
  X <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, d)
  th <- fit$opt$par
  s2 <- exp(2 * th[["betad"]])
  lt <- 2 * (th[["theta"]] - th[["betad"]])
  gl <- reml_gls(X, d$y, d$g, exp(lt))
  # J: the estimate moves with log sd and against log sigma alike
  h <- 1e-4
  db <- (reml_gls(X, d$y, d$g, exp(lt + 2 * h))$beta -
           reml_gls(X, d$y, d$g, exp(lt - 2 * h))$beta) / (2 * h)
  J <- cbind(-db, db)
  Vt <- allow_warnings(vcov(fit, full = TRUE),
                       "a product such as a' V a")
  Vt <- Vt[nrow(Vt) - 1:0, nrow(Vt) - 1:0]
  # the variance of a' beta: conditional part through the factor, plus
  # what the variance parameters' uncertainty adds
  vref <- function(A) {
    s2 * rowSums((A %*% gl$Ri)^2) + rowSums((A %*% J %*% Vt) * (A %*% J))
  }
  V <- allow_warnings(vcov(fit), "a product such as a' V a")
  se <- unname(sqrt(diag(V)))
  expect_lt(max(abs(se / sqrt(vref(diag(ncol(X)))) - 1)), 1e-4)
  nd <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  A <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, nd)
  fp <- fitted(fit, newdata = nd, re_formula = NA)
  expect_lt(max(abs(unname(fp[, "Est.Error"]) / sqrt(vref(A)) - 1)), 1e-4)
  # the estimates are the GLS ones at the REML variance parameters
  expect_lt(max(abs(fixef(fit)[, "Estimate"] - gl$beta) / se), 1e-3)
})

test_that("a well-conditioned REML fit keeps its natural coordinates", {
  # the guard's absent case: degree 2 is far below the condition number
  # at which the design's coordinates are taken
  d <- reml_data(1, 2)
  fit <- frm(y ~ x1 + x2 + (1 | g), data = d, REML = TRUE)
  expect_null(fit$beta_whiten)
  expect_null(sdr_of(fit)$joint_root)
  expect_null(frmtmb:::beta_whitener(fit$frame))
})

# the restricted log-likelihood at (log sigma, log sd) = th, in double
# through the whitened QR: -2 l = (n - p) log(2 pi) + log|V| +
# log|X' V^-1 X| + r' V^-1 r
reml_loglik <- function(X, y, g, th) {
  s2 <- exp(2 * th[1])
  t <- exp(2 * (th[2] - th[1]))
  Xt <- X
  yt <- y
  ld <- 0
  for (ix in split(seq_along(y), g)) {
    ng <- length(ix)
    a <- (1 - 1 / sqrt(1 + ng * t)) / ng
    yt[ix] <- y[ix] - a * sum(y[ix])
    Xt[ix, ] <- X[ix, , drop = FALSE] -
      a * rep(colSums(X[ix, , drop = FALSE]), each = ng)
    ld <- ld + log(1 + ng * t)
  }
  q <- qr(Xt, tol = 1e-12)
  n <- length(y)
  p <- ncol(X)
  rss <- sum(qr.resid(q, yt)^2)
  -((n - p) * log(2 * pi) + n * log(s2) + ld +
      2 * sum(log(abs(diag(qr.R(q))))) - p * log(s2) + rss / s2) / 2
}

test_that("REML logLik of a whitened fit is the restricted likelihood", {
  # Laplace over gamma integrated over beta / |det T|: 0.71.0's round-1
  # build of this lane put logLik 13.5, 20.3 and 28.4 below lme4's at
  # degrees 4 to 6 (review of lane illse2, finding 1). Against the
  # restricted log-likelihood in double (whitened QR) at the fit's own
  # estimates, and lme4's where lme4 keeps its digits (degree 4 and 5;
  # at degree 6 lme4's own criterion is 1.6e-4 off, dev/illse2-fx-remlll.R)
  skip_if_not_installed("lme4")
  for (deg in 4:6) {
    d <- reml_data(1, deg)
    xn <- paste0("x", seq_len(deg))
    fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                  "+ (1 | g)"))
    fit <- frm(fo, data = d, REML = TRUE)
    expect_false(is.null(fit$beta_whiten))
    X <- stats::model.matrix(stats::as.formula(paste("~",
                                       paste(xn, collapse = " + "))), d)
    ll <- as.numeric(stats::logLik(fit))
    ref <- reml_loglik(X, d$y, d$g, unname(fit$opt$par))
    expect_lt(abs(ll / ref - 1), 1e-9)
    if (deg <= 5) {
      l <- suppressWarnings(lme4::lmer(fo, data = d, REML = TRUE))
      expect_lt(abs(ll / as.numeric(stats::logLik(l)) - 1), 1e-7)
    }
  }
})

test_that("a shifted covariate leaves the REML logLik where it was", {
  # the restricted likelihood does not depend on a shift of a covariate;
  # `Days + 1e5` takes the whitener and `Days` does not, and they gave
  # -587.808 and -586.742 for one model
  skip_if_not_installed("lme4")
  d <- lme4::sleepstudy
  d$s <- d$Days + 1e5
  f0 <- frm(Reaction ~ Days + (1 | Subject), data = d, REML = TRUE)
  f1 <- frm(Reaction ~ s + (1 | Subject), data = d, REML = TRUE)
  expect_null(f0$beta_whiten)
  expect_false(is.null(f1$beta_whiten))
  l0 <- as.numeric(stats::logLik(f0))
  expect_lt(abs(as.numeric(stats::logLik(f1)) / l0 - 1), 1e-9)
})
