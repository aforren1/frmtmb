# Multi-dimensional gp() terms (exact + Hilbert-space) and kriging
# prediction for the exact form.

gp2_data <- function() {
  set.seed(23)
  n <- 80
  grid <- expand.grid(x1 = seq(0, 4, by = 0.5), x2 = seq(0, 4, by = 0.5))
  pu <- as.matrix(grid[sample(nrow(grid), 50), ])
  D1 <- outer(pu[, 1], pu[, 1], "-")
  D2 <- outer(pu[, 2], pu[, 2], "-")
  Ku <- 1.2^2 * exp(-D1^2 / (2 * 1.5^2) - D2^2 / (2 * 0.8^2))
  uu <- drop(crossprod(chol(Ku + diag(1e-8, 50)), rnorm(50)))
  idx <- sample(50, n, replace = TRUE)
  data.frame(y = 1 + uu[idx] + rnorm(n, 0, 0.4),
             x1 = pu[idx, 1], x2 = pu[idx, 2])
}

# unique coordinate rows in the frame's (lexicographic) order, plus the
# pieces the closed-form references need
gp2_ref_parts <- function(dg) {
  posdf <- unique(data.frame(x1 = dg$x1, x2 = dg$x2))
  posdf <- posdf[order(posdf$x1, posdf$x2), ]
  P <- as.matrix(posdf)
  list(
    P = P,
    Dm1 = outer(P[, 1], P[, 1], "-")^2,
    Dm2 = outer(P[, 2], P[, 2], "-")^2,
    Zi = outer(paste(dg$x1, dg$x2), paste(P[, 1], P[, 2]), `==`) * 1
  )
}

test_that("2-D exact gp() matches direct GP marginal ML", {
  dg <- gp2_data()
  n <- nrow(dg)
  rp <- gp2_ref_parts(dg)
  np <- nrow(rp$P)

  f2 <- frm(bf(y ~ gp(x1, x2)) + gaussian(), data = dg)
  # anisotropic (brms default): one lengthscale per dimension
  nll <- function(p) {
    Kp <- exp(2 * p[2]) * (exp(-rp$Dm1 / (2 * exp(2 * p[3])) -
                                rp$Dm2 / (2 * exp(2 * p[4]))) +
                             diag(1e-6, np))
    S <- rp$Zi %*% Kp %*% t(rp$Zi) + exp(2 * p[5]) * diag(n)
    r <- dg$y - p[1]
    0.5 * (determinant(S)$modulus + sum(r * solve(S, r)) +
             n * log(2 * pi))
  }
  op <- stats::optim(c(1, log(1.2), log(1.5), log(0.8), log(0.4)), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(f2)) + op$value), 1e-4)
  cv <- confint_varcorr(f2)
  expect_setequal(cv$term,
                  c("sd(gp)", "range(gp, x1)", "range(gp, x2)"))

  # iso = TRUE: one shared lengthscale over Euclidean distance
  fi <- frm(bf(y ~ gp(x1, x2, iso = TRUE)) + gaussian(), data = dg)
  nlli <- function(p) {
    Kp <- exp(2 * p[2]) * (exp(-(rp$Dm1 + rp$Dm2) / (2 * exp(2 * p[3]))) +
                             diag(1e-6, np))
    S <- rp$Zi %*% Kp %*% t(rp$Zi) + exp(2 * p[4]) * diag(n)
    r <- dg$y - p[1]
    0.5 * (determinant(S)$modulus + sum(r * solve(S, r)) +
             n * log(2 * pi))
  }
  opi <- stats::optim(c(1, log(1.2), log(1), log(0.4)), nlli,
                      method = "BFGS",
                      control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fi)) + opi$value), 1e-4)
  expect_setequal(confint_varcorr(fi)$term, c("sd(gp)", "range(gp)"))
  # iso nests inside anisotropic
  expect_gte(as.numeric(logLik(f2)), as.numeric(logLik(fi)) - 1e-6)

  # guards: dimension cap, basis-size cap
  expect_error(frm(bf(y ~ gp(x1, x2, x1, x2)) + gaussian(), data = dg),
               "at most 3 dimensions")
  expect_error(frm(bf(y ~ gp(x1, x2, k = 40)) + gaussian(), data = dg),
               "cap 1000")
})

test_that("exact gp() kriging matches the closed form", {
  set.seed(17)
  n <- 120
  xg <- round(runif(n, 0, 10), 1)
  pos <- sort(unique(xg))
  Dm <- abs(outer(pos, pos, "-"))
  K <- 1.5^2 * exp(-Dm^2 / 8)
  u <- drop(crossprod(chol(K + diag(1e-8, length(pos))),
                      rnorm(length(pos))))
  dg <- data.frame(y = 2 + u[match(xg, pos)] + rnorm(n, 0, 0.5), x = xg)
  fg <- frm(bf(y ~ gp(x)) + gaussian(), data = dg)

  # observed positions: exact reproduction of fitted values and the
  # pre-kriging se.fit (indicator fast path)
  p_obs <- predict(fg, newdata = dg[1:8, ], se.fit = TRUE)
  expect_equal(unname(p_obs$fit), unname(fitted(fg)[1:8]),
               tolerance = 1e-12)

  # fitted kernel pieces (nugget included, as in the fit)
  bk <- fg$frame$re_blocks[[1]]
  est <- fg$estimates
  th <- est$theta[bk$theta_idx]
  sd2 <- exp(2 * th[1])
  rho <- exp(th[2])
  Kfit <- sd2 * (exp(-Dm^2 / (2 * rho^2)) + diag(1e-6, length(pos)))
  bhat <- est$b[bk$b_idx]

  # grid of unseen positions (data is on a 0.1 grid; offsets miss it)
  xs <- seq(0.05, 9.55, by = 0.5)
  expect_false(any(xs %in% pos))
  nd <- data.frame(x = xs)
  pk <- predict(fg, newdata = nd, se.fit = TRUE)

  # conditional mean K* K^-1 b_hat plus the fixed-effect part
  Kst <- sd2 * exp(-outer(xs, pos, "-")^2 / (2 * rho^2))
  Xr <- t(solve(Kfit, t(Kst)))
  mu_hand <- est$beta[1] + drop(Xr %*% bhat)
  expect_lt(max(abs(pk$fit - mu_hand)), 1e-8)

  # variance decomposition: delta-method part over (beta, b) plus the
  # GP conditional variance diag(K** - K* K^-1 K*')
  jc <- frmtmb:::get_joint_cov(fg)
  cp <- c(which(jc$names == "beta")[1], which(jc$names == "b")[bk$b_idx])
  A <- cbind(1, Xr)
  V <- jc$V[cp, cp]
  v_delta <- rowSums((A %*% V) * A)
  ev_hand <- sd2 * (1 + 1e-6) - rowSums(Xr * Kst)
  expect_lt(max(abs(pk$se.fit^2 - (v_delta + ev_hand))), 1e-8)

  # the GP conditional variance is strictly positive at every unseen
  # position, so the kriging se strictly exceeds the pure-interpolation
  # (delta-method-only) se there. It need NOT exceed the se at the
  # nearest observed position: interpolation smooths the b uncertainty
  # (Var(E[f*|b]) shrinks between knots) faster than the conditional
  # variance grows near a knot.
  expect_true(all(ev_hand > 0))
  expect_true(all(pk$se.fit > sqrt(v_delta)))

  # mixed newdata (observed + unseen rows) goes through the kriging
  # path; observed rows still reproduce the indicator behavior
  ndm <- data.frame(x = c(dg$x[1:3], 0.05))
  pm <- predict(fg, newdata = ndm, se.fit = TRUE)
  expect_lt(max(abs(pm$fit[1:3] - p_obs$fit[1:3])), 1e-8)
  expect_lt(max(abs(pm$se.fit[1:3] - p_obs$se.fit[1:3])), 1e-8)
  expect_lt(abs(pm$fit[4] - pk$fit[1]), 1e-8)
  expect_lt(abs(pm$se.fit[4] - pk$se.fit[1]), 1e-8)
})

test_that("gp() k/c/iso resolve in the formula environment", {
  dg <- gp2_data()
  kk <- 10
  cc <- 1.5

  f_lit <- frm(bf(y ~ gp(x1, x2, k = 10, c = 1.5)) + gaussian(), data = dg)
  f_var <- frm(bf(y ~ gp(x1, x2, k = kk, c = cc)) + gaussian(), data = dg)
  expect_equal(as.numeric(logLik(f_var)), as.numeric(logLik(f_lit)))
  expect_equal(coef(f_var), coef(f_lit))
  # the spec keeps the evaluated scalars, so newdata prediction needs no
  # access to kk/cc
  ndg <- expand.grid(x1 = c(0.25, 2.25), x2 = c(0.75, 2.75))
  expect_equal(predict(f_var, newdata = ndg),
               predict(f_lit, newdata = ndg))

  # arbitrary expressions, not just names
  f_expr <- frm(bf(y ~ gp(x1, x2, k = 5 + 5, c = 1.5)) + gaussian(),
                data = dg)
  expect_equal(as.numeric(logLik(f_expr)), as.numeric(logLik(f_lit)))

  use_iso <- TRUE
  f_iso <- frm(bf(y ~ gp(x1, x2, k = 4, iso = use_iso)) + gaussian(),
               data = dg)
  expect_true(f_iso$frame$re_blocks[[1]]$gp_iso)
  expect_equal(as.numeric(logLik(f_iso)),
               as.numeric(logLik(frm(bf(y ~ gp(x1, x2, k = 4,
                                               iso = TRUE)) + gaussian(),
                                     data = dg))))

  expect_error(frm(bf(y ~ gp(x1, x2, k = no_such_k)) + gaussian(),
                   data = dg),
               "cannot evaluate k")
  expect_error(frm(bf(y ~ gp(x1, x2, k = c(4, 5))) + gaussian(),
                   data = dg),
               "must be a single value")
  expect_error(frm(bf(y ~ gp(x1, x2, c = -1)) + gaussian(), data = dg),
               "must be positive")
  expect_error(frm(bf(y ~ gp(x1, x2, iso = 1)) + gaussian(), data = dg),
               "must be TRUE or FALSE")
})

test_that("2-D Hilbert-space gp() approximates the exact fit", {
  dg <- gp2_data()
  f2 <- frm(bf(y ~ gp(x1, x2)) + gaussian(), data = dg)
  # brms's convention rescales both coordinates by one shared factor (the
  # largest pairwise distance over the whole coordinate matrix), so the
  # boundary is L = c per dimension. The wider effective domain lets the
  # default c = 1.25 reach 0.3 logLik at k = 10, where the pre-brms
  # half-range convention needed c = 1.5 for less accuracy.
  fh <- frm(bf(y ~ gp(x1, x2, k = 10)) + gaussian(), data = dg)
  expect_equal(fh$frame$linpreds[["y.mu"]]$gps[[1]]$L, c(1.25, 1.25),
               tolerance = 1e-12)
  expect_lt(abs(as.numeric(logLik(fh)) - as.numeric(logLik(f2))), 0.3)

  # unseen grid (positions are multiples of 0.5): hsgp basis evaluates
  # anywhere, the exact fit kriges; the surfaces agree
  ndg <- expand.grid(x1 = c(0.25, 1.25, 2.25, 3.25),
                     x2 = c(0.75, 1.75, 2.75))
  p_e <- predict(f2, newdata = ndg)
  p_h <- predict(fh, newdata = ndg, se.fit = TRUE)
  expect_lt(max(abs(p_h$fit - p_e)), 0.06)
  expect_true(all(is.finite(p_h$se.fit)))

  # in-sample newdata reproduces the fit exactly (stored scaling)
  expect_equal(unname(predict(fh, newdata = dg)), unname(fitted(fh)),
               tolerance = 1e-12)

  # the lengthscales are estimated on the rescaled inputs but reported in
  # data units, so they track the exact fit's per-dimension ranges
  cvh <- confint_varcorr(fh)
  cve <- confint_varcorr(f2)
  expect_setequal(cvh$term,
                  c("sd(gp)", "range(gp, x1)", "range(gp, x2)"))
  for (tm in c("range(gp, x1)", "range(gp, x2)")) {
    expect_lt(abs(cvh$estimate[cvh$term == tm] /
                    cve$estimate[cve$term == tm] - 1), 0.15)
  }

  # c = is per covariate in brms; a vector widens one boundary only
  fv <- frm(bf(y ~ gp(x1, x2, k = 8, c = c(1.5, 2))) + gaussian(),
            data = dg, dry_run = "frame")
  expect_equal(fv$linpreds[["y.mu"]]$gps[[1]]$L, c(1.5, 2),
               tolerance = 1e-12)
  expect_error(frm(bf(y ~ gp(x1, x2, k = 8, c = c(1, 2, 3))) + gaussian(),
                   data = dg, dry_run = "frame"),
               "length 1 or the number of variables")
})

test_that("rr(d=) and se(sigma=) resolve in the formula environment", {
  set.seed(31)
  n <- 120
  g <- factor(rep(1:12, each = 10))
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- rnorm(n, 1 + 0.5 * x1, 1)
  dd <- data.frame(y = y, x1 = x1, x2 = x2, g = g)
  rk <- 1L
  f_var <- suppressWarnings(
    frm(bf(y ~ x1 + rr(x1 + x2 | g, d = rk)) + gaussian(), data = dd))
  f_lit <- suppressWarnings(
    frm(bf(y ~ x1 + rr(x1 + x2 | g, d = 1)) + gaussian(), data = dd))
  expect_equal(logLik(f_var), logLik(f_lit), tolerance = 1e-8)
  expect_error(
    frm(bf(y ~ x1 + rr(x1 + x2 | g, d = no_rank)) + gaussian(),
        data = dd),
    "cannot evaluate d")
  expect_error(
    frm(bf(y ~ x1 + rr(x1 + x2 | g, d = 0)) + gaussian(), data = dd),
    "positive whole number")

  dse <- data.frame(yi = rnorm(20, 0, 0.6), sei = runif(20, 0.3, 0.5))
  flag <- TRUE
  fs_var <- frm(bf(yi | se(sei, sigma = flag) ~ 1) + gaussian(),
                data = dse)
  fs_lit <- frm(bf(yi | se(sei, sigma = TRUE) ~ 1) + gaussian(),
                data = dse)
  expect_equal(logLik(fs_var), logLik(fs_lit), tolerance = 1e-8)
  expect_error(
    frm(bf(yi | se(sei, sigma = maybe) ~ 1) + gaussian(), data = dse),
    "cannot evaluate sigma")
})
