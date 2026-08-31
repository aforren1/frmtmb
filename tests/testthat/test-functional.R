# Function-on-scalar regression in mixed-model form: curves stacked long,
# coefficient functions as smooths (varying-coefficient s(t, by = x)),
# curve-level random intercepts. The pffr/mgcv representation.
test_that("function-on-scalar regression matches mgcv", {
  set.seed(141)
  n_sub <- 40
  nt <- 20
  tt <- seq(0, 1, length.out = nt)
  xs <- rnorm(n_sub)
  b_sub <- rnorm(n_sub, 0, 0.4)
  f0 <- function(t) sin(2 * pi * t)
  f1 <- function(t) cos(pi * t)
  dd <- expand.grid(t = tt, id = seq_len(n_sub))
  dd$x <- xs[dd$id]
  dd$y <- f0(dd$t) + dd$x * f1(dd$t) + b_sub[dd$id] +
    rnorm(nrow(dd), 0, 0.25)
  dd$id <- factor(dd$id)

  fit <- frm(bf(y ~ s(t) + s(t, by = x) + (1 | id)) + gaussian(),
             data = dd)
  ref <- mgcv::gam(y ~ s(t) + s(t, by = x) + s(id, bs = "re"),
                   data = dd, method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - (-as.numeric(ref$gcv.ubre))),
            1e-3)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 0.05)

  # the estimated coefficient function beta1(t) = d eta / d x at fixed t:
  # predict at x = 1 minus x = 0, population level
  nd1 <- data.frame(t = tt, x = 1, id = factor(1, levels = levels(dd$id)))
  nd0 <- data.frame(t = tt, x = 0, id = factor(1, levels = levels(dd$id)))
  beta1_hat <- predict(fit, newdata = nd1, re.form = NA) -
    predict(fit, newdata = nd0, re.form = NA)
  expect_lt(max(abs(beta1_hat - f1(tt))), 0.25)
})
