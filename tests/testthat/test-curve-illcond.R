## A curve on an ill-conditioned fixed design: frmtmb (>= the lane
## illse core) holds that fit's covariance as a factor, and the grid
## covariance is formed through it, since C V C' from the whole matrix
## loses every digit there.

sp_poly_fit <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  fo <- y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) +
    I(x^8)
  list(d = d, fo = fo, fit = frmtmb::frm(fo, data = d))
}

test_that("a curve on a degree-8 raw polynomial has lm()'s bands", {
  o <- sp_poly_fit()
  g <- data.frame(x = seq(1, 2, length.out = 11))
  cv <- frm_curve(o$fit, newdata = g, simultaneous = FALSE)
  lmf <- stats::lm(o$fo, data = o$d)
  n <- nrow(o$d)
  k <- length(stats::coef(lmf))
  ref <- stats::predict(lmf, newdata = g, se.fit = TRUE)$se.fit *
    sqrt((n - k) / n)
  expect_lt(max(abs(cv$.se / ref - 1)), 1e-4)
  # the matrix formed from V whole: the check against frm_linpred()
  # refused it, or its diagonal was not the standard error
  expect_lt(max(abs(sqrt(diag(attr(cv, "Sigma"))) / cv$.se - 1)), 1e-6)
  # a contrast between two grids, against the same difference by lm()'s
  # QR
  g2 <- data.frame(x = g$x + 0.05)
  dif <- frm_curve(o$fit, newdata = g, contrast = g2, simultaneous = FALSE)
  X1 <- stats::model.matrix(stats::delete.response(stats::terms(lmf)), g)
  X2 <- stats::model.matrix(stats::delete.response(stats::terms(lmf)), g2)
  Ri <- backsolve(qr.R(lmf$qr), diag(k))
  s <- sqrt(sum(stats::residuals(lmf)^2) / n)
  dref <- sqrt(rowSums(((X1 - X2) %*% Ri)^2)) * s
  expect_lt(max(abs(dif$.se / dref - 1)), 1e-4)
  # the first derivative by its central difference of design rows, to
  # the stencil's own accuracy (dev/illse-splderiv.R: within 4e-4)
  gd <- data.frame(x = seq(1.1, 1.9, length.out = 9))
  d1 <- frm_curve_deriv(o$fit, var = "x", order = 1, newdata = gd,
                        simultaneous = FALSE)
  D <- outer(gd$x, 0:8, function(x, j) j * x^pmax(j - 1, 0))
  ddref <- sqrt(rowSums((D %*% Ri)^2)) * s
  expect_lt(max(abs(d1$.se / ddref - 1)), 1e-3)
})
