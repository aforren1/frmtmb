# Tracer validation for quadseam-short-hetrace.R: the lane's own case
# (poissd2 seed 29, quadrature, default grad_tol), where the fix round
# says the Gauss-Kronrod tape's AD Hessian is NaN. Not a package test.
test_that("tracer sees a non-finite he()", {
  set.seed(29)
  d <- data.frame(g = factor(rep(1:50, each = 3)), x = rnorm(150))
  u <- rnorm(50, 0, 2)
  d$y <- rpois(nrow(d), exp(-1 + 0.3 * d$x + u[d$g]))
  f <- suppressWarnings(frm(y ~ x + (1 | g), family = poisson(), data = d,
                            quadrature = TRUE))
  H <- f$obj$he(f$opt$par)
  cat("\nVALIDATE he finite:", all(is.finite(H)), "\n")
  v <- frmtmb:::grad_verdict(f, control = list(grad_tol = 1e-12))
  cat("VALIDATE forced headroom:", v$headroom, "\n")
  expect_true(TRUE)
})
