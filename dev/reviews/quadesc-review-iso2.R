# Review of lane quadesc: narrowing the NaN of group 11 at x*
# (quadesc-review-iso-ref.txt). 1. RTMB pnorm(log.p = TRUE) derivative on a
# fine log grid of x in [-1e14, -1e2]. 2. Group 11's log integrand
# l(u; theta) taped in theta at fixed u over a wide u grid: where its
# gradient is not finite, and which observation's term carries it.
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
tp <- RTMB::MakeTape(function(x) RTMB::pnorm(x, log.p = TRUE), 0)
J <- tp$jacfun()
xs <- -10^seq(2, 14, length.out = 4001)
dj <- vapply(xs, function(x) J(x)[1], 1)
bad <- !is.finite(dj)
cat(sprintf("pnorm log deriv: %d of %d non-finite on [-1e14, -1e2]\n",
            sum(bad), length(xs)))
if (any(bad)) {
  cat("  non-finite at x in [", format(min(xs[bad])), ",",
      format(max(xs[bad])), "]\n")
}
rel <- abs(dj / -xs - 1)
for (p in c(4, 5, 6, 7, 8, 9, 10, 12, 14)) {
  i <- which.min(abs(xs + 10^p))
  cat(sprintf("  x = -1e%d deriv %.6g rel err vs -x %.3g\n", p, dj[i],
              rel[i]))
}

sn <- function(n, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  dl * (abs(stats::rnorm(n)) - sqrt(2 / pi)) +
    sqrt(1 - dl^2) * stats::rnorm(n)
}
set.seed(1011)
invisible(stats::rnorm(60 + 15 + 60))
invisible(stats::rnorm(50 + 10))
invisible(sn(50, 1.5))
d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:16, each = 5)))
d$y <- 0.5 * d$x + stats::rnorm(16, 0, 1)[d$g] + sn(80, 2)
d$y <- d$y * (1 + 2e-15 * 9)
xstar <- c(-4.2234454e-01, 3.9856199e-01, -4.3553262e-01, 1.9925412e+06,
           -5.7424988e+06, -3.4749739e-01)
dg <- d[d$g == "11", ]
cat("group 11 x:", format(dg$x, digits = 4), "\n")
al <- xstar[4] + xstar[5] * dg$x
cat("group 11 alpha_i at x*:", format(al, digits = 4), "\n")
for (obs in c(0, seq_len(nrow(dg)))) {
  f <- function(th) {
    u <- th[7]
    mu <- th[1] + th[2] * dg$x + u
    a <- th[4] + th[5] * dg$x
    l <- RTMBdist::dskewnorm2(dg$y, mu, exp(th[3]), a, log = TRUE)
    if (obs > 0) l <- l[obs]
    sum(l) + RTMB::dnorm(u, 0, exp(th[6]), log = TRUE)
  }
  tpg <- RTMB::MakeTape(f, c(xstar, 0))
  Jg <- tpg$jacfun()
  us <- c(-rev(10^seq(-3, 4, length.out = 2000)),
          10^seq(-3, 4, length.out = 2000))
  nf <- vapply(us, function(u) {
    jj <- Jg(c(xstar, u))
    !all(is.finite(jj[1:6]))
  }, NA)
  cat(sprintf("obs %s: gradient non-finite at %d of %d u values%s\n",
              if (obs == 0) "all" else obs, sum(nf), length(us),
              if (any(nf)) sprintf(", u in [%.4g, %.4g]",
                                   min(us[nf]), max(us[nf])) else ""))
  if (obs > 0 && any(nf)) {
    u <- us[nf][1]
    dl <- al[obs] / sqrt(1 + al[obs]^2)
    om <- exp(xstar[3]) / sqrt(1 - 2 * dl^2 / pi)
    xi <- xstar[1] + xstar[2] * dg$x[obs] + u - om * dl * sqrt(2 / pi)
    z <- (dg$y[obs] - xi) / om
    cat(sprintf("   first such u %.6g: z %.6g alpha*z %.6g value %.6g\n",
                u, z, al[obs] * z, tpg(c(xstar, u))))
  }
}
cat("END\n")
