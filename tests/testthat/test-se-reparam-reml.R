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

test_that("a whitened REML fit with a lost sd predicts through the factor", {
  skip_on_cran()
  # Backlog item 3370. With no group effect the sd runs to 0 and is lost,
  # so the joint covariance is repaired (joint_cov_repair()), which
  # inverted the inner block in natural coordinates, where beta's block
  # has the design's condition number squared. vcov() kept its digits;
  # fitted(re_formula = NA) was 0.7 to 3.2 percent off at degree 7
  # (seeds 1 to 5, dev/selimit2-3370.R). With the sd at 0 the GLS weights
  # are 1, so the conditional covariance is s2 (X'X)^-1, and the variance
  # parameters move the estimate by nothing
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  y0 <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  X <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, d)
  nd <- data.frame(x = c(1, 1.25, 1.5, 1.75, 2), g = d$g[1])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  A <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, nd)
  lost <- 0L
  for (k in 0:19) {
    d$y <- y0 * (1 + 2e-15 * k)
    fit <- frm(fo, data = d, REML = TRUE)
    expect_false(is.null(fit$beta_whiten))
    lost <- lost + identical(unname(sdr_of(fit)$se_lost), "boundary")
    th <- fit$opt$par
    gl <- reml_gls(X, d$y, d$g, exp(2 * (th[["theta"]] - th[["betad"]])))
    vref <- exp(2 * th[["betad"]]) * rowSums((A %*% gl$Ri)^2)
    fp <- fitted(fit, newdata = nd, re_formula = NA)
    expect_lt(max(abs(unname(fp[, "Est.Error"]) / sqrt(vref) - 1)), 1e-6)
  }
  expect_gt(lost, 0L)
})

test_that("a whitened REML fit above the size limit reads sparse solves", {
  skip_on_cran()
  # Backlog item 4952. Above options(frmtmb.joint_root_max) no dense
  # factor is formed. Through 0.82.0 the natural inverse then gave
  # standard errors of exactly 0 or 7 percent high at 4510 joint
  # parameters, with no warning, and 0.83.0 made them NaN with a
  # warning. They now come from sparse solves in the design's
  # coordinates. The limit is lowered before the fit, which forms its
  # factor at its first standard error, so the fit stays small (30 joint
  # parameters). The sd runs to 0 and is lost, so the joint covariance
  # is repaired and the conditional covariance is the exact one
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  y0 <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  X <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, d)
  nd <- data.frame(x = c(1, 1.25, 1.5, 1.75, 2), g = d$g[1])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  A <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, nd)
  withr::local_options(frmtmb.joint_root_max = 20)
  lost <- 0L
  for (k in 0:19) {
    d$y <- y0 * (1 + 2e-15 * k)
    fit <- frm(fo, data = d, REML = TRUE)
    expect_false(is.null(fit$beta_whiten))
    lost <- lost + identical(unname(sdr_of(fit)$se_lost), "boundary")
    expect_true(inherits(get_joint_cov(fit)$V, "frm_joint_sparse"))
    expect_identical(attr(sdr_of(fit)$joint_root, "rows"),
                     which(rownames(sdr_of(fit)$jointPrecision) != "b"))
    th <- fit$opt$par
    gl <- reml_gls(X, d$y, d$g, exp(2 * (th[["theta"]] - th[["betad"]])))
    s2 <- exp(2 * th[["betad"]])
    vref <- s2 * rowSums((A %*% gl$Ri)^2)
    fp <- fitted(fit, newdata = nd, re_formula = NA)
    expect_lt(max(abs(unname(fp[, "Est.Error"]) / sqrt(vref) - 1)), 1e-6)
    V <- allow_warnings(vcov(fit), "a product such as a' V a")
    expect_lt(max(abs(sqrt(unname(diag(V))) /
                        sqrt(s2 * rowSums(gl$Ri^2)) - 1)), 1e-6)
  }
  expect_gt(lost, 0L)
})

test_that("the sparse route agrees with the dense one on every reader", {
  # One fit read twice: through the dense factor (the default limit) and
  # through sparse solves (the limit lowered before the fit). Both are the
  # same algebra, so they agree to the rounding of the factors: within
  # 5.6e-10 at 3510 joint parameters (dev/sparsese-agree.R). The group
  # effect is real, so nothing is lost and the whole inverse is read
  d <- reml_data(1, 7)
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  nd <- data.frame(x = seq(1, 2, by = 0.25), g = d$g[1:5])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  lossy <- "a product such as a' V a"
  read <- function(lim) {
    withr::local_options(frmtmb.joint_root_max = lim)
    fit <- frm(fo, data = d, REML = TRUE)
    list(fit = fit,
         noRE = fitted(fit, newdata = nd, re_formula = NA)[, "Est.Error"],
         RE = fitted(fit, newdata = nd)[, "Est.Error"],
         rows = fitted(fit)[, "Est.Error"],
         lin = frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit,
         V = allow_warnings(vcov(fit), lossy),
         Vf = allow_warnings(vcov(fit, full = TRUE), lossy),
         jc = allow_warnings(frm_joint_cov(fit)$V, lossy))
  }
  dn <- read(joint_root_max)
  sp <- read(20)
  expect_null(attr(sdr_of(dn$fit)$joint_root, "rows"))
  expect_false(is.null(attr(sdr_of(sp$fit)$joint_root, "rows")))
  expect_identical(dn$fit$opt$par, sp$fit$opt$par)
  expect_length(sdr_of(sp$fit)$se_lost, 0L)
  for (nm in c("noRE", "RE", "rows", "lin")) {
    expect_lt(max(abs(unname(sp[[nm]]) / unname(dn[[nm]]) - 1)), 1e-7,
              label = nm)
  }
  # matrices, against their largest entry
  for (nm in c("V", "Vf", "jc")) {
    a <- unname(as.matrix(sp[[nm]]))
    b <- unname(as.matrix(dn[[nm]]))
    expect_lt(max(abs(a - b)) / max(abs(b)), 1e-9, label = nm)
  }
})

test_that("the item's whole-route fit above 4000 parameters is exact", {
  skip_on_cran()
  # Backlog item 4952, seed 1: 4500 groups of 2 rows, 4510 joint
  # parameters, the default limit. The sd is not lost, so nothing is
  # repaired; vcov() and fitted() read the whole inverse, which 0.82.0
  # and 0.83.0 took in the natural coordinates: fitted(re_formula = NA)
  # at x = 1 read 0.023272941 against a conditional lower bound of
  # 0.024359899, 4.5 percent low, and vcov() the same, with no warning
  set.seed(1)
  ng <- 4500L
  n <- 2L * ng
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(seq_len(ng), length.out = n)),
                  f = factor(rep(c("a", "b"), length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  fit <- frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d,
             REML = TRUE)
  expect_gt(nrow(sdr_of(fit)$jointPrecision), joint_root_max)
  X <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, d)
  nd <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  A <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, nd)
  th <- fit$opt$par
  s2 <- exp(2 * th[["betad"]])
  lt <- 2 * (th[["theta"]] - th[["betad"]])
  gl <- reml_gls(X, d$y, d$g, exp(lt))
  vref <- s2 * rowSums((A %*% gl$Ri)^2)
  lossy <- "a product such as a' V a"
  if (!length(sdr_of(fit)$se_lost)) {
    # what the variance parameters' uncertainty adds, as in the first
    # test of this file; a rounding that loses the sd leaves the
    # conditional part alone
    h <- 1e-4
    db <- (reml_gls(X, d$y, d$g, exp(lt + 2 * h))$beta -
             reml_gls(X, d$y, d$g, exp(lt - 2 * h))$beta) / (2 * h)
    J <- cbind(-db, db)
    Vt <- allow_warnings(vcov(fit, full = TRUE), lossy)
    Vt <- Vt[nrow(Vt) - 1:0, nrow(Vt) - 1:0]
    vref <- vref + rowSums((A %*% J %*% Vt) * (A %*% J))
  }
  fp <- fitted(fit, newdata = nd, re_formula = NA)
  expect_lt(max(abs(unname(fp[, "Est.Error"]) / sqrt(vref) - 1)), 1e-4)
  # vcov() is formed from the sparse whole inverse (joint_cov_whole());
  # the rows of the factor sdr_unwhiten() keeps, which hypothesis(), the
  # draws and interop read, give the same standard errors and matrix
  root <- sdr_of(fit)$joint_root
  expect_false(is.null(attr(root, "rows")))
  jb <- which(rownames(sdr_of(fit)$jointPrecision) == "beta")
  expect_lt(max(abs(sqrt(rowSums((A %*% root[jb, ])^2)) / sqrt(vref) - 1)),
            1e-4)
  V <- allow_warnings(vcov(fit), lossy)
  expect_lt(max(abs(unname(V) - tcrossprod(root[jb, ]))) / max(abs(V)),
            1e-12)
})

# 30 groups of 10 rows, a group effect, a raw polynomial: a gaussian
# response and a three-category ordinal one
sparse_route_data <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(300, 1, 2),
                  g = factor(rep(1:30, length.out = 300)),
                  f = factor(rep(c("a", "b"), length.out = 300)))
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  eta <- sin(3 * d$x) + stats::rnorm(30, 0, 0.5)[d$g] + 0.3 * (d$f == "b")
  d$y <- eta + stats::rnorm(300, 0, 0.3)
  d$yc <- factor(cut(0.8 * eta + stats::rlogis(300), c(-Inf, 0, 0.8, Inf),
                     labels = c("lo", "mid", "hi")), ordered = TRUE)
  d
}

test_that("frm_curve() and frm_lp_basis()'s V_root read the sparse route", {
  # Review of lane sparsese, F1. frm_lp_basis() gave no V_root on the
  # sparse route, so frm_curve() refused ("disagrees with
  # frm_linpred(se.fit = TRUE) by 0.00132 relative") and the documented
  # recipe through V_root was not there. With the factor taken in the
  # whitened coordinates both agree with the dense route to 7.6e-10
  # (dev/sparsese-fix1.R, seeds 1 to 3)
  d <- sparse_route_data(1)
  fo <- y ~ f + x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  grid <- data.frame(x = seq(1, 2, by = 0.1), f = factor("a", c("a", "b")))
  for (k in 1:7) grid[[paste0("x", k)]] <- grid$x^k
  read <- function(lim) {
    withr::local_options(frmtmb.joint_root_max = lim)
    fit <- frm(fo, data = d, REML = TRUE)
    lb <- frm_lp_basis(fit, newdata = grid, re_formula = NA)
    list(fit = fit, lb = lb,
         se = frm_linpred(fit, newdata = grid, re_formula = NA,
                          se.fit = TRUE)$se.fit,
         cv = frm_curve(fit, newdata = grid, re_formula = NA,
                        simultaneous = FALSE)$.se)
  }
  dn <- read(joint_root_max)
  sp <- read(20)
  expect_true(inherits(get_joint_cov(sp$fit)$V, "frm_joint_sparse"))
  expect_false(is.null(sp$lb$V_root))
  rec <- sqrt(rowSums((sp$lb$A %*% sp$lb$V_root)^2) +
                rowSums((sp$lb$A %*% sp$lb$V_rest) * sp$lb$A))
  expect_lt(max(abs(rec / sp$se - 1)), 1e-7)
  expect_lt(max(abs(sp$cv / dn$cv - 1)), 1e-7)
  expect_lt(max(abs(sp$se / dn$se - 1)), 1e-7)
})

test_that("a cumulative fit's fitted() reads the sparse route", {
  # Review of lane sparsese, F2. fit_fd_se() needs the factor's rows at
  # the group effects (fd_joint_root()); the sparse route had none, and
  # the per-parameter differences that took over lose every digit of
  # J V J' on a raw polynomial: Est.Error up to 7.39 for a probability,
  # 124 times the dense route's. The factor is now taken over the drawn
  # parameters in the whitened coordinates (jcs_fd_root()): within 4e-7
  # of the dense route (dev/sparsese-fix1.R, seeds 1 to 3)
  d <- sparse_route_data(1)
  nd <- d[c(1:6, 151:156), ]
  read <- function(lim) {
    withr::local_options(frmtmb.joint_root_max = lim)
    fit <- frm(yc ~ x1 + x2 + x3 + x4 + x5 + (1 | g), data = d,
               family = cumulative(), REML = TRUE)
    list(fit = fit, new = fitted(fit, newdata = nd)[, "Est.Error", ],
         rows = fitted(fit)[, "Est.Error", ])
  }
  dn <- read(joint_root_max)
  sp <- read(20)
  expect_false(is.null(sp$fit$beta_whiten))
  expect_true(inherits(get_joint_cov(sp$fit)$V, "frm_joint_sparse"))
  expect_identical(dn$fit$opt$par, sp$fit$opt$par)
  for (nm in c("new", "rows")) {
    expect_lt(max(abs(sp[[nm]] / dn[[nm]] - 1)), 1e-5, label = nm)
  }
  expect_length(sdr_of(sp$fit)$se_lost, 0L)
  expect_lt(max(sp$rows), 0.5)
})

test_that("a whitened REML fit whose sparse route fails keeps NaN SEs", {
  # The sparse route fails only where the inner block of the whitened
  # precision does not factor; the natural inverse is then no answer,
  # and the coefficients' standard errors stay NaN with a warning, as in
  # 0.83.0. Mocked here: the block of a fit at a maximum factors
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  withr::local_options(frmtmb.joint_root_max = 20)
  fit <- frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d,
             REML = TRUE)
  expect_identical(unname(sdr_of(fit)$se_lost), "boundary")
  local_mocked_bindings(jcs_joint_cov = function(...) NULL)
  nd <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  fp <- allow_warnings(fitted(fit, newdata = nd, re_formula = NA),
                       "does not factor",
                       require = "above 20 joint parameters")
  expect_true(all(is.nan(fp[, "Est.Error"])))
  expect_true(all(is.finite(fp[, "Estimate"])))
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
