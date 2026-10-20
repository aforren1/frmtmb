# Routes that invert the covariance of an ill-conditioned fixed design:
# cooks.distance() and the Wald rules of anova() on frm_multiple() fits
# read the factor of the covariance (root_inv_quad()), where solve() of
# the whole matrix stops as computationally singular. Lane illse2,
# dev/illse2-findings.md, item 1.

inv_conds <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

# a degree-8 raw polynomial on [1, 2], the design of test-se-reparam.R
inv_poly_data <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  d
}
inv_fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
inv_xn <- paste0("x", 1:8)

test_that("cooks.distance() reads the factor on a degree-8 design", {
  # The reference is the exact gaussian Hessian H in (beta, log sigma)
  # at the fit's own point: V^-1 = H, and d' H d formed through X d
  # keeps about kappa(X) epsilon (within 3e-9 of a 256-bit reference,
  # dev/illse2-p1ref.R). frmtmb 0.71.0 stopped with "system is
  # computationally singular" (dev/illse2-p1.R)
  d <- inv_poly_data(1)
  fit <- frm(inv_fo, data = d)
  expect_false(is.null(sdr_of(fit)$cov_fixed_root))
  inf <- influence(fit, groups = "g", data = d)
  cd <- inv_conds(cooks.distance(inf))
  expect_length(cd$warnings, 0L)
  X <- stats::model.matrix(inv_fo, d)
  p <- fit$opt$par
  r <- d$y - drop(X %*% p[1:9])
  s2 <- exp(2 * p[10])
  D <- sweep(inf$fixed, 2, inf$fixed_full)
  ref <- apply(D, 1L, function(di) {
    u <- drop(X %*% di[1:9])
    (sum(u^2) + 4 * di[10] * sum(u * r) + 2 * di[10]^2 * sum(r^2)) / s2
  }) / ncol(D)
  expect_lt(max(abs(cd$value / ref - 1)), 1e-4)
})

test_that("the Wald rules of anova() read the factors over every power", {
  # Three imputations that differ in y only, so every fit has the same
  # design. With the cross term of beta and log sigma at the maximum
  # neglected (1e-10 here, dev/illse2-p1ref.R), U_j over x1..x8 is
  # s2_j (Xc' Xc)^-1 with Xc the centered columns, and every quadratic
  # form of D1 and D2 is a sum of squares of Xc times a vector
  d <- inv_poly_data(2)
  dl <- lapply(1:3, function(j) {
    set.seed(2000 + j)
    dj <- d
    dj$y <- dj$y + stats::rnorm(nrow(dj), 0, 0.05)
    dj
  })
  mm <- frm_multiple(inv_fo, data = dl)
  Xc <- scale(stats::model.matrix(inv_fo, d)[, -1L], scale = FALSE)
  Q <- vapply(mm$fits, function(f) f$opt$par[2:9], numeric(8))
  s2 <- vapply(mm$fits, function(f) exp(2 * f$opt$par[[10]]), 0)
  form <- function(v, s) sum(drop(Xc %*% v)^2) / s
  qb <- rowMeans(Q)
  tr <- sum(apply((Q - qb) / sqrt(2), 2L, form, s = mean(s2)))
  r1 <- max(0, (1 + 1 / 3) * tr / 8)
  st1 <- form(qb, mean(s2)) / ((1 + r1) * 8)
  dj <- vapply(1:3, function(j) form(Q[, j], s2[j]), 0)
  r2 <- max(0, (1 + 1 / 3) * stats::var(sqrt(dj)))
  st2 <- (mean(dj) / 8 - 2 * r2) / (1 + r2)
  a1 <- inv_conds(anova(mm, method = "D1", constraint = inv_xn))
  expect_length(a1$warnings, 0L)
  expect_lt(abs(a1$value$statistic / st1 - 1), 1e-4)
  expect_lt(abs(a1$value$riv / r1 - 1), 1e-4)
  a2 <- inv_conds(anova(mm, method = "D2", use = "wald",
                        constraint = inv_xn))
  expect_length(a2$warnings, 0L)
  expect_lt(abs(a2$value$statistic / st2 - 1), 1e-4)
})

test_that("a factor whose rows are dependent gives no inverse", {
  # the case the routes turn into NaN with a warning, never into a
  # number solve() would make of a singular matrix
  K <- rbind(c(1, 2, 3), c(2, 4, 6))
  expect_null(frmtmb:::root_inv_quad(K, c(1, 1)))
  # the absent case: independent rows give the form of K K'
  K2 <- rbind(c(1, 2, 3), c(0, 1, 5))
  v <- c(0.3, -1.2)
  expect_equal(frmtmb:::root_inv_quad(K2, v),
               drop(crossprod(v, solve(tcrossprod(K2), v))))
})

test_that("a well-conditioned design keeps the inverse it had", {
  # no factor, so cooks.distance() and D1 form solve(V) exactly as
  # before: the guard's absent case, bit for bit
  set.seed(3)
  d <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(60)
  fit <- frm(y ~ x, data = d)
  expect_null(sdr_of(fit)$cov_fixed_root)
  expect_null(frmtmb:::coef_vcov_root(fit))
  inf <- influence(fit, groups = "g", data = d)
  V <- frmtmb:::interop_vcov(fit)
  Vi <- solve(V)
  D <- sweep(inf$fixed, 2, inf$fixed_full)
  ref <- vapply(seq_len(nrow(D)), function(i) {
    drop(t(D[i, ]) %*% Vi %*% D[i, ]) / ncol(V)
  }, 0)
  expect_identical(unname(cooks.distance(inf)), ref)
})
