# bernoulli b1 * exp(b2 * x) from the default start: stationary, not a saddle,
# below brms's best (dev/nlfix-saddle.R, seed 9205, the same draws). Which
# point is it?
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                   c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
set.seed(9205)
d$yg <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
d$yc <- rpois(n, 1.5 * exp(0.6 * d$x))
d$yzi <- ifelse(runif(n) < 0.25, 0L, rpois(n, 1.5 * exp(0.6 * d$x)))
d$yhu <- ifelse(runif(n) < 0.25, 0L, 1L + rpois(n, 1.5 * exp(0.6 * d$x)))
d$ypo <- rgamma(n, 3, 3 / (1.5 * exp(0.6 * d$x)))
d$ybi <- rbinom(n, 1, plogis(1.5 * exp(0.6 * d$x) - 2))
f <- bf(ybi ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
a <- suppressWarnings(frm(f, family = bernoulli(), data = d))
cat("default", as.numeric(logLik(a)), "conv", a$opt$convergence, "beta",
    a$estimates$beta, "maxgrad", max(abs(a$obj$gr(a$opt$par))), "\n")
for (s in list(c(.3, .3), c(-.3, .3), c(.3, -.3), c(-.3, -.3), c(-1, .5), c(1, .5))) {
  r <- suppressWarnings(frm(f, family = bernoulli(), data = d, start = list(
    beta = c("b1_(Intercept)" = s[1], "b2_(Intercept)" = s[2]))))
  cat("start", s, as.numeric(logLik(r)), "conv", r$opt$convergence, "beta",
      r$estimates$beta, "\n")
}
