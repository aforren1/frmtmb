# Extracted from test-v16.R:149

# test -------------------------------------------------------------------------
set.seed(73)
n <- 150
z <- rnorm(n)
x <- rnorm(n, 0.5 + 0.8 * z, 0.7)
y <- rnorm(n, 1 + 0.6 * x + 0.3 * z, 0.9)
d2 <- data.frame(y = y, x = x, z = z)
f2 <- frm(bf(y ~ mi(x) + z) + gaussian() +
              bf(x | mi() ~ z) + gaussian(), data = d2)
fy <- frm(bf(y ~ x + z) + gaussian(), data = d2)
fx <- frm(bf(x ~ z) + gaussian(), data = d2)
expect_lt(abs(as.numeric(logLik(f2)) -
                  as.numeric(logLik(fy)) - as.numeric(logLik(fx))),
            1e-6)
d3 <- d2
d3$x[1:10] <- NA
d3$z[1:5] <- NA
f3 <- frm(bf(y ~ mi(x) + z) + gaussian() +
              bf(x | mi() ~ z) + gaussian(), data = d3)
expect_equal(nobs(f3), n - 5L)
expect_length(f3$estimates$miss, 5L)
expect_error(frm(bf(y ~ mi(x) + z) + gaussian(), data = d2),
               "matching imputation model")
expect_error(frm(bf(y ~ mi(x) + z) + gaussian() +
                     bf(x | mi() ~ z) + poisson(), data = d3),
               "gaussian or student")
expect_error(frm(bf(y ~ mi(x) * z) + gaussian() +
                     bf(x | mi() ~ z) + gaussian(), data = d2),
               "standalone")
