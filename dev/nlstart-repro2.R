# The help's bernoulli example of a second maximum (dev/nlfix-bern.R, seeds
# 9202 and 9205), and the cost of the probe on bodies with three and four
# coefficients it varies (prod3, two_exp), in one process: min of 3
# interleaved rounds, off2 / off the control, each time a block of 30 fits
# over 30.
# Usage: NLSTART_ARM=before|after Rscript dev/nlstart-repro2.R
source("dev/nlstart-common.R")
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
cat("bernoulli help example, default:", ll(tt(frm(f, family = bernoulli(),
                                                 data = d))), "\n")
cat("  start (1, 0.5):", ll(tt(frm(f, family = bernoulli(), data = d,
  start = list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.5))))),
  "\n")
if (arm == "before") quit("no")
# three coefficients the probe varies: a logistic curve (SSlogis's form) and
# two exponentials
set.seed(7304)
n <- 300
d3 <- data.frame(x = runif(n, 0, 10), z = rnorm(n))
d3$y <- 5 / (1 + exp((4 - d3$x) / 1.5)) + rnorm(n, 0, 0.3)
d3$y2 <- 1.2 * exp(0.3 * d3$x / 5) + 0.8 * exp(-0.5 * d3$z) +
  rnorm(n, 0, 0.3)
designs <- list(
  prod3 = list(f = bf(y ~ b1 * exp(b2 * x / 5) * (1 + b3 * z),
                      b1 + b2 + b3 ~ 1, nl = TRUE), d = d3),
  two_exp = list(f = bf(y2 ~ b1 * exp(b2 * x / 5) + b3 * exp(b4 * z),
                        b1 + b2 + b3 + b4 ~ 1, nl = TRUE), d = d3))
fit_t <- function(ds, probe) {
  old <- options(frmtmb.nl_sign_probe = probe)
  on.exit(options(old))
  t0 <- proc.time()[["elapsed"]]
  for (k in 1:30) r <- tt(frm(ds$f, data = ds$d))
  list(t = (proc.time()[["elapsed"]] - t0) / 30, r = r)
}
for (nm in names(designs)) {
  ds <- designs[[nm]]
  arms <- c(off = FALSE, on = TRUE, off2 = FALSE)
  tm <- list(off = numeric(0), on = numeric(0), off2 = numeric(0))
  last <- list()
  for (round in 0:2) {
    ord <- names(arms)[(seq_along(arms) + round) %% 3L + 1L]
    for (a in ord) {
      z <- fit_t(ds, arms[[a]])
      tm[[a]] <- c(tm[[a]], z$t)
      last[[a]] <- z$r
    }
  }
  m <- vapply(tm, min, 0)
  cat(sprintf("K3 %s min s: off %.2f on %.2f off2 %.2f | on/off %.3f control %.3f\n",
              nm, m[["off"]], m[["on"]], m[["off2"]], m[["on"]] / m[["off"]],
              m[["off2"]] / m[["off"]]))
  cat("   off:", ll(last$off), "\n   on: ", ll(last$on), "\n")
}
