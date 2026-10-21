# Reviewer, fix round 1: a poisson fit truncated to [31, 60], where the
# normalizer's window (30, 60] has its lower edge at the last count of
# the summed form, against an independent ML fit with stats::ppois().
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
set.seed(8)
n <- 4000
x <- rnorm(n)
y <- rpois(n, exp(3.3 + 0.4 * x))
keep <- y >= 31 & y <= 60
d <- data.frame(x = x[keep], y = y[keep])
cat("rows", nrow(d), "\n")
nll <- function(p) {
  mu <- exp(p[1] + p[2] * d$x)
  -sum(dpois(d$y, mu, log = TRUE) - log(ppois(60, mu) - ppois(30, mu)))
}
o <- nlminb(c(3.3, 0.4), nll)
cat(sprintf("independent: %s logLik %.6f\n", paste(signif(o$par, 8), collapse = " "), -o$objective))
fit <- tryCatch(frm(y | trunc(lb = 31, ub = 60) ~ x, family = poisson(), data = d),
                error = function(e) conditionMessage(e))
if (is.character(fit)) cat("frmtmb ERROR", fit, "\n") else {
  cat(sprintf("frmtmb:      %s logLik %.6f conv %d\n", paste(signif(fit$opt$par, 8), collapse = " "),
              as.numeric(logLik(fit)), fit$opt$convergence))
  cat(sprintf("frmtmb objective at the independent optimum %.6f vs %.6f\n", -fit$obj$fn(o$par), -o$objective))
}
cat("== rows whose normalizer is wrong at the independent optimum\n")
mu <- exp(o$par[1] + o$par[2] * d$x)
ns <- asNamespace("frmtmb")
got <- ns$row_lpdf(ns$as_frmtmb_family(poisson()), d$y, d$y, list(mu = mu),
                   list(trunc_lb = 31, trunc_ub = 60), NULL)
ref <- dpois(d$y, mu, log = TRUE) - log(ppois(60, mu) - ppois(30, mu))
cat("max |row - ref|", max(abs(got - ref)), " rows off by > 1e-6:", sum(abs(got - ref) > 1e-6), "of", nrow(d), "\n")
cat("== the objective along the intercept (slope 0.4): frmtmb vs independent\n")
for (b0 in c(3.3, 2.5, 2, 1.5, 1.2, 1, 0.5)) {
  cat(sprintf("intercept %.1f: frmtmb -fn %.4f  independent %.4f\n", b0, -fit$obj$fn(c(b0, 0.4)), -nll(c(b0, 0.4))))
}
cat("== the same with the window's mass as a log-sum of dpois(31:60) (exact)\n")
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
truell <- function(p) {
  mu <- exp(p[1] + p[2] * d$x)
  lw <- vapply(mu, function(m) lse(dpois(31:60, m, log = TRUE)), 0)
  sum(dpois(d$y, mu, log = TRUE) - lw)
}
for (b0 in c(2, 1.5, 1, 0.5)) cat(sprintf("intercept %.1f: frmtmb -fn %.4f  exact %.4f\n", b0, -fit$obj$fn(c(b0, 0.4)), truell(c(b0, 0.4))))
