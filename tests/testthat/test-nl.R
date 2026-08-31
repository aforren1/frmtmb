test_that("nonlinear fixed-effects model matches nls", {
  set.seed(111)
  n <- 200
  x <- runif(n, 0, 5)
  y <- 2.5 * exp(-0.7 * x) + rnorm(n, 0, 0.15)
  dd <- data.frame(y = y, x = x)

  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) + gaussian(),
             data = dd, start = list(beta = c(1, 0.3)))
  ref <- nls(y ~ a * exp(-b * x), data = dd, start = list(a = 1, b = 0.3))

  expect_lt(abs(fixef(fit)$a[[1]] - coef(ref)[["a"]]), 1e-4)
  expect_lt(abs(fixef(fit)$b[[1]] - coef(ref)[["b"]]), 1e-4)
  # ML sigma^2 = RSS/n at the same coefficients
  sig_ml <- sqrt(sum(residuals(ref)^2) / n)
  expect_lt(abs(exp(fixef(fit)$sigma[[1]]) - sig_ml), 1e-4)
})

test_that("nonlinear mixed model matches a hand-rolled reference", {
  set.seed(112)
  n_g <- 25; n_per <- 20
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- runif(n_g * n_per, 0, 5)
  a_g <- 2.5 + rnorm(n_g, 0, 0.5)
  y <- a_g[g] * exp(-0.7 * x) + rnorm(length(x), 0, 0.15)
  dd <- data.frame(y = y, x = x, g = g)

  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) +
               gaussian(),
             data = dd, start = list(beta = c(2, 0.5)))

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  nll_ref <- function(p) {
    nll <- -sum(RTMB::dnorm(p$u, 0, exp(p$lsd), log = TRUE))
    a <- p$a0 + p$u[gi]
    mu <- a * exp(-p$b0 * xv)
    nll - sum(RTMB::dnorm(yv, mu, exp(p$ls), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(a0 = 2, b0 = 0.5, ls = 0, lsd = 0,
                              u = numeric(n_g)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)

  # nlpar random effects show up in ranef and VarCorr
  expect_length(VarCorr(fit), 1)
  expect_identical(dim(ranef(fit)[[1]]), c(as.integer(n_g), 1L))
})

test_that("nl prediction and post-processing", {
  set.seed(113)
  n <- 150
  x <- runif(n, 0, 5)
  dd <- data.frame(y = 2 * exp(-0.5 * x) + rnorm(n, 0, 0.1), x = x)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) + gaussian(),
             data = dd, start = list(beta = c(1, 0.3)))

  expect_equal(predict(fit, newdata = dd), predict(fit), tolerance = 1e-8)
  expect_equal(fitted(fit), predict(fit, type = "response"),
               tolerance = 1e-8)
  a_hat <- predict(fit, dpar = "a")
  expect_lt(stats::sd(a_hat), 1e-10)   # intercept-only nlpar is constant
  nd <- data.frame(x = c(0, 1, 2))
  p <- predict(fit, newdata = nd)
  expect_equal(p[1], fixef(fit)$a[[1]], tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_error(predict(fit, se.fit = TRUE), "se.fit is not supported")
  expect_length(residuals(fit), n)
})

test_that("nl validation errors are clear", {
  expect_error(bf(y ~ a * exp(-b * x), nl = TRUE), "parameter formula")
  expect_error(frm(bf(y ~ a * exp(-b * x), a ~ 1, cc ~ 1, nl = TRUE) +
                     gaussian(),
                   data = NULL, dry_run = "spec"),
               "not used in the model formula")
})
