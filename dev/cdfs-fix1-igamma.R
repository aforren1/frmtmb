# Fix round 1, finding 1: log_igamma_two() across the whole shape range
# after the exponent cap, against stats::pgamma() in both tails, and the
# two fits the review found failing, against an independent ML fit on
# stats::pgamma() / stats::ppois().
#   Rscript dev/cdfs-fix1-igamma.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
cat("== accuracy: relative error of the log, floored at one\n")
shapes <- c(0.01, 0.1, 0.5, 1, 2, 5, 12, 30, 60, 80, 100, 300, 1e3, 3e3,
            1e4, 1.25e4, 1.3e4, 3e4, 1e5, 3e5, 1e6, 1e7)
worst_all <- 0
nan_all <- 0
for (s in shapes) {
  sd <- sqrt(s)
  x <- unique(sort(c(s * c(1e-3, 0.01, 0.1, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8,
                           0.85, 0.9, 0.95, 0.99, 1, 1.01, 1.05, 1.5, 2, 4),
                     pmax(s + sd * c(-40, -20, -10, -5, -2, -1, 0, 1, 2, 5,
                                     10, 20, 40), 1e-3 * s))))
  x <- x[x > 0]
  p <- ns$log_igamma_two(x, s)
  rP <- stats::pgamma(x, s, log.p = TRUE)
  rQ <- stats::pgamma(x, s, lower.tail = FALSE, log.p = TRUE)
  ok <- is.finite(rP) & is.finite(rQ)
  e <- pmax(abs(p$lP - rP) / pmax(1, abs(rP)),
            abs(p$lQ - rQ) / pmax(1, abs(rQ)))[ok]
  nn <- sum(is.nan(p$lP) | is.nan(p$lQ))
  worst_all <- max(worst_all, e)
  nan_all <- nan_all + nn
  cat(sprintf("a %-8g points %3d NaN %d worst %.2e\n", s, sum(ok), nn, max(e)))
}
cat(sprintf("ALL: NaN %d worst %.2e\n", nan_all, worst_all))

cat("== gradient finite at a = 1e5, x = 0.7 a and 0.3 a\n")
f <- function(p) {
  g <- ns$log_igamma_two(exp(p[2]), exp(p[1]))
  g$lP + g$lQ
}
for (r in c(0.3, 0.7)) {
  tp <- RTMB::MakeTape(f, c(log(1e5), log(r * 1e5)))
  cat(sprintf("x/a %.1f: %s\n", r,
              paste(signif(tp$jacobian(c(log(1e5), log(r * 1e5))), 8),
                    collapse = " ")))
}

cat("== censored Gamma at shape 1e5 (review's data)\n")
set.seed(3)
d <- data.frame(x = rnorm(200))
d$y <- rgamma(200, shape = 1e5, rate = 1e5 / exp(1 + 0.1 * d$x))
qq <- quantile(d$y, 0.8)
d$cc <- as.integer(d$y > qq)
d$y[d$cc == 1] <- qq
fit <- frm(y | cens(cc) ~ x, family = Gamma(link = "log"), data = d)
nll_g <- function(p) {
  mu <- exp(p[1] + p[2] * d$x)
  sh <- exp(p[3])
  -sum(ifelse(d$cc == 1,
              stats::pgamma(d$y, sh, sh / mu, lower.tail = FALSE, log.p = TRUE),
              stats::dgamma(d$y, sh, sh / mu, log = TRUE)))
}
ind <- stats::nlminb(c(1, 0.1, log(1e5)), nll_g)
cat(sprintf("frmtmb logLik %.8f conv %d | independent %.8f | diff %.2e\n",
            as.numeric(logLik(fit)), fit$opt$convergence, -ind$objective,
            as.numeric(logLik(fit)) + ind$objective))

cat("== poisson right censored at 75000 (counts near 1e5), brms's P(Y > k)\n")
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
d$y <- rpois(n, exp(11 + 0.3 * d$x))
d$cc <- 0
d$cc[d$y > 75000] <- 1
d$y[d$cc == 1] <- 75000
cat("censored rows:", sum(d$cc), "\n")
fit <- frm(y | cens(cc) ~ x, family = poisson(), data = d)
nll_p <- function(p) {
  mu <- exp(p[1] + p[2] * d$x)
  -sum(ifelse(d$cc == 1,
              stats::ppois(d$y, mu, lower.tail = FALSE, log.p = TRUE),
              stats::dpois(d$y, mu, log = TRUE)))
}
ind <- stats::nlminb(c(11, 0.3), nll_p)
cat(sprintf("frmtmb logLik %.6f conv %d | independent %.6f | diff %.2e\n",
            as.numeric(logLik(fit)), fit$opt$convergence, -ind$objective,
            as.numeric(logLik(fit)) + ind$objective))
cat("q 100000 mu 70000 lS:",
    ns$poisson_cdf_log(1e5, 7e4)$lS, "ref",
    stats::ppois(1e5, 7e4, lower.tail = FALSE, log.p = TRUE), "\n")
cat("END\n")
