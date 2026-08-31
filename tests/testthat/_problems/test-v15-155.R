# Extracted from test-v15.R:155

# test -------------------------------------------------------------------------
set.seed(6)
n <- 400
inc <- sample(0:3, n, replace = TRUE)
x <- rnorm(n)
y <- 1 + c(0, 0.5, 0.8, 1)[inc + 1] * 2 + 0.3 * x + rnorm(n, 0, 0.8)
dd <- data.frame(y = y, inc = inc, x = x)
fit <- frm(bf(y ~ x + mo(inc)) + gaussian(), data = dd)
nll <- function(p) {
    zr <- exp(c(0, p[4:6]))
    cz0 <- c(0, cumsum(zr / sum(zr)))
    -sum(stats::dnorm(y, p[1] + p[2] * x + p[3] * 3 * cz0[inc + 1],
                      exp(p[7]), log = TRUE))
  }
op <- stats::optim(c(1, 0.3, 0.6, 0, 0, 0, log(0.8)), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)
p_new <- predict(fit, newdata = data.frame(x = 0, inc = 0:3))
expect_true(all(diff(p_new) >= -1e-8))
expect_equal(unname(fitted(fit)),
               unname(predict(fit, newdata = dd, type = "response")),
               tolerance = 1e-10)
ps <- predict(fit, newdata = data.frame(x = 0, inc = 0:3),
                se.fit = TRUE)
expect_true(all(is.finite(ps$se.fit)))
dd$incf <- factor(inc, levels = 0:3, ordered = TRUE)
ff <- frm(bf(y ~ x + mo(incf)) + gaussian(), data = dd)
expect_loglik_equal(ff, fit, tol = 1e-6)
expect_error(frm(bf(y ~ mo(inc) * x) + gaussian(), data = dd),
               "standalone")
