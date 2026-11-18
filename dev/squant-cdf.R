# Lane squant, item 4(e): the accuracy of quantile_elf()'s log CDF pair
# against stats::pbeta(log.p = TRUE), over shapes from the grid of
# smoothing x quantile, and on the tape against the doubles.
#
# Run: Rscript dev/squant-cdf.R > dev/squant-out/cdf.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")

ts <- c(-680, -300, -100, -40, -10, -3, -1, -0.1, 0, 0.1, 1, 3, 10, 40,
        100, 300, 680)
worst <- 0
cat("h\tp\ta\tb\tmax_rel_lF\tmax_rel_lS\tmax_complement\n")
for (h in c(0.02, 0.1, 0.25, 0.5, 1, 4, 20)) {
  for (p in c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) {
    a <- h * (1 - p)
    b <- h * p
    q <- ts * h  # mu 0, sigma 1: t = q / h
    pr <- ns$elf_cdf_log(q, 0, 1, p, 1 - p, h)
    # each tail from the side where its argument is below one half, so
    # that plogis() never rounds to one
    lo <- stats::plogis(-abs(ts))
    rf <- ifelse(ts <= 0, stats::pbeta(lo, a, b, log.p = TRUE),
                 stats::pbeta(lo, b, a, lower.tail = FALSE, log.p = TRUE))
    rs <- ifelse(ts >= 0, stats::pbeta(lo, b, a, log.p = TRUE),
                 stats::pbeta(lo, a, b, lower.tail = FALSE, log.p = TRUE))
    # relative error, floored at one in the log value
    ef <- abs(pr$lF - rf) / pmax(1, abs(rf))
    es <- abs(pr$lS - rs) / pmax(1, abs(rs))
    ec <- abs(exp(pr$lF) + exp(pr$lS) - 1)
    cat(sprintf("%.2f\t%.2f\t%.4f\t%.4f\t%.2e\t%.2e\t%.2e\n", h, p, a, b,
                max(ef), max(es), max(ec)))
    worst <- max(worst, ef, es)
  }
}
cat("WORST relative error in the log CDF pair:", format(worst, digits = 3),
    "\n")

# beyond |t| = 690: the continued tails against the leading term of the
# incomplete beta, a t - log(a B(a, b)) and -b t - log(b B(a, b))
h <- 0.25
p <- 0.9
a <- h * (1 - p)
b <- h * p
tt <- c(-5000, -1000, -700)
far <- ns$elf_cdf_log(tt * h, 0, 1, p, 1 - p, h)
cat("far lF:", format(far$lF, digits = 17), "\n")
cat("lead  :", format(a * tt - log(a) - lbeta(a, b), digits = 17), "\n")
tt <- c(700, 1000, 5000)
far <- ns$elf_cdf_log(tt * h, 0, 1, p, 1 - p, h)
cat("far lS:", format(far$lS, digits = 17), "\n")
cat("lead  :", format(-b * tt - log(b) - lbeta(a, b), digits = 17), "\n")

# the tape against the doubles, value and gradient in mu
f <- function(par) {
  v <- ns$elf_cdf_log(c(-3, 0.2, 5), par[1], exp(par[2]), 0.9, 0.1, 0.25)
  sum(v$lF) + sum(v$lS)
}
tp <- RTMB::MakeTape(f, c(0.1, 0))
cat("tape", format(tp(c(0.1, 0)), digits = 17), "doubles",
    format(f(c(0.1, 0)), digits = 17), "\n")
g <- tp$jacfun()(c(0.1, 0))
e <- 1e-6
gn <- c((f(c(0.1 + e, 0)) - f(c(0.1 - e, 0))) / (2 * e),
        (f(c(0.1, e)) - f(c(0.1, -e))) / (2 * e))
cat("gradient tape", format(g, digits = 10), "finite difference",
    format(gn, digits = 10), "\n")
