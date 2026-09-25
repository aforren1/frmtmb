# Reviewer 2, item 7: the three "tape defects" in plain RTMB, no frmtmb.
suppressPackageStartupMessages(library(RTMB))
options(digits = 10, width = 150)
cat("RTMB", as.character(packageVersion("RTMB")), " TMB",
    as.character(packageVersion("TMB")), "\n")

cat("\n(1) d/dx pnorm(x, log.p = TRUE). Truth phi(x)/Phi(x), about -x for x << 0\n")
tp <- MakeTape(function(x) pnorm(x, log.p = TRUE), 0)
tpu <- MakeTape(function(x) pnorm(x, lower.tail = FALSE, log.p = TRUE), 0)
truth <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
for (x in c(-10, -1e3, -3e4, -1e5, -3e5, -1e6, -2e6, -1e7, -2e8, -8.3e9, 38, 1e5, 8.3e9)) {
  g <- tp$jacobian(x)
  cat(sprintf("x = %10.3g  value %14.8g (R %14.8g)  tape d %14.8g  truth %14.8g  rel %9.2e\n",
              x, tp(x), pnorm(x, log.p = TRUE), g, truth(x), abs(g - truth(x)) / abs(truth(x))))
}
cat("upper-tail spelling, pnorm(x, lower.tail = FALSE, log.p = TRUE), truth -phi(x)/(1-Phi(x))\n")
truthu <- function(x) -exp(dnorm(x, log = TRUE) - pnorm(x, lower.tail = FALSE, log.p = TRUE))
for (x in c(1e3, 1e5, 1e6, 2e8)) {
  g <- tpu$jacobian(x)
  cat(sprintf("x = %10.3g  tape d %14.8g  truth %14.8g  rel %9.2e\n", x, g, truthu(x),
              abs(g - truthu(x)) / abs(truthu(x))))
}
cat("second derivative of pnorm(x, log.p = TRUE): tape vs truth (about -1 for x << 0)\n")
h <- tp$jacfun()
for (x in c(-10, -1e3, -1e5, -1e6)) cat(sprintf("x = %9.3g  tape d2 %14.8g\n", x, h$jacobian(x)))

cat("\n(2) zero drift: d2/dx2 of sqrt(x^2 + 1e-20) at 0. Exact value 1/sqrt(1e-20) = 1e10\n")
t2 <- MakeTape(function(x) sqrt(x * x + 1e-20), 0)
cat("tape d2 at 0:", t2$jacfun()$jacobian(0), "  (this is the regularizer's true curvature, ",
    "not a tape error)\n")
t3 <- MakeTape(function(x) log(sinh(x) / x), 1)
cat("log(sinh(x)/x) tape d2 at x = 1e-3:", t3$jacfun()$jacobian(1e-3), " truth ~ 1/3\n")

cat("\n(3) tanh: tape d and d2 against truth (d = sech^2, d2 = -2 tanh sech^2)\n")
t4 <- MakeTape(function(x) tanh(x), 0)
h4 <- t4$jacfun()
for (x in c(20, 40, 300, 700, 710, 711, 800, 1e4, -711)) {
  cat(sprintf("x = %7g  d %12.5g  d2 %12.5g\n", x, t4$jacobian(x), h4$jacobian(x)))
}
cat("same for cosh overflow proxy: d2 of log(cosh(x)) at 711:",
    MakeTape(function(x) log(cosh(x)), 0)$jacfun()$jacobian(711), "\n")
