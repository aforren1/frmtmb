# check_convergence() at a reported optimum whose gradient is not
# finite. The gradient test compares a finite number with grad_tol, so
# such a fit used to pass in silence: lane ordmix's review found a
# three-component ordinal mixture that nlminb reported converged
# (X-convergence) with an infinite gradient and no warning
# (dev/ordmix-rev-probit-sat.R).

test_that("a non-finite gradient at the optimum warns, naming it", {
  set.seed(20261070)
  d <- data.frame(x = rnorm(60))
  d$y <- 1 + 0.5 * d$x + rnorm(60)
  fit <- frm(y ~ x, data = d, family = gaussian())
  # the fit as it is: a finite gradient, nothing to say
  expect_no_warning(frmtmb:::check_convergence(fit, fit$control))
  # the same fit, its gradient made infinite in the slope
  bad <- fit
  g0 <- fit$obj$gr
  bad$obj$gr <- function(p) {
    g <- g0(p)
    g[2] <- Inf
    g
  }
  expect_warning(frmtmb:::check_convergence(bad, fit$control),
                 "The gradient at the reported optimum is not finite (beta)",
                 fixed = TRUE)
  bad$obj$gr <- function(p) {
    g <- g0(p)
    g[1] <- NaN
    g
  }
  expect_warning(frmtmb:::check_convergence(bad, fit$control),
                 "is not finite", fixed = TRUE)
})

test_that("the three-component mixture the review found now warns", {
  # the data of dev/ordmix-rev-probit-sat.R: the unused draws keep the
  # random stream where the review had it
  set.seed(20261005 + 31)
  n <- 400
  x <- rnorm(n)
  rnorm(n)
  sample(c("a", "b"), n, TRUE)
  rbinom(n, 1, 0.4)
  rlogis(n)
  runif(n)
  runif(n, 0.5, 2)
  cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
  lat <- c(1.5, -0.8, 0.3)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
  d <- data.frame(y = as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf))),
                  x = x)
  fit <- allow_warnings(
    frm(bf(y ~ x), family = mixture(cumulative("probit"), sratio("cloglog"),
                                    acat()), data = d),
    c("gradient at the reported optimum is not finite",
      "degenerate boundary"),
    require = "gradient at the reported optimum is not finite")
  expect_false(all(is.finite(fit$obj$gr(fit$opt$par))))
})
