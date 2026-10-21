# Reviewer: isolate the asym_laplace truncated-gradient mismatch:
# upper bound only, lower only, both; constant vs column; and the sum
# of log F(ub) over rows taped by hand with a vector mu.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(77)
n <- 200
x <- rnorm(n)
y <- rnorm(n)
keep <- y >= -1.5 & y <= 2
d <- data.frame(x = x[keep], y = y[keep], lbv = -1.5, ubv = 2)
chk <- function(fo) {
  ob <- frm(fo, family = asym_laplace(), data = d, dry_run = "objective")$obj
  p <- ob$par
  g <- as.numeric(ob$gr(p))
  fd <- numDeriv::grad(ob$fn, p)
  cat(sprintf("%-46s AD %s | FD %s\n", deparse(fo$formula), paste(signif(g, 6), collapse = ","),
              paste(signif(fd, 6), collapse = ",")))
}
chk(bf(y | trunc(ub = ubv) ~ x, sigma ~ 1, quantile ~ 1))
chk(bf(y | trunc(lb = lbv) ~ x, sigma ~ 1, quantile ~ 1))
chk(bf(y | trunc(ub = 2) ~ x, sigma ~ 1, quantile ~ 1))
chk(bf(y | trunc(ub = 2) ~ 1, sigma ~ 1, quantile ~ 1))
f <- ns$as_frmtmb_family(asym_laplace())
cat("== by hand: sum log F(2) over rows, mu = b0 + b1 x\n")
h <- function(p) {
  mu <- p[1] + p[2] * d$x
  sum(log(f$lcdf(rep(2, nrow(d)), list(mu = mu, sigma = exp(p[3]), quantile = 1 / (1 + exp(-p[4]))), list())))
}
p <- c(0.05, 0.02, 0.1, 0.2)
cat("AD", signif(RTMB::MakeTape(h, p)$jacobian(p), 7), "\nFD", signif(numDeriv::grad(h, p), 7), "\n")
h2 <- function(p) {
  mu <- p[1] + p[2] * d$x
  u <- (2 - mu) / exp(p[3])
  sum(ns$sign_pos(u))
}
cat("sum sign_pos(u) AD", signif(RTMB::MakeTape(h2, p)$jacobian(p), 7), "\n")
h3 <- function(p) {
  mu <- p[1] + p[2] * d$x
  u <- (2 - mu) / exp(p[3])
  sum(ns$ad_min(u, 0)) + sum(ns$ad_max(u, 0))
}
cat("sum ad_min + ad_max AD", signif(RTMB::MakeTape(h3, p)$jacobian(p), 7),
    " FD", signif(numDeriv::grad(h3, p), 7), "\n")
cat("== the same objective one step off the start (mu intercept + 1e-3)\n")
ob <- frm(bf(y | trunc(ub = 2) ~ 1, sigma ~ 1, quantile ~ 1), family = asym_laplace(), data = d,
          dry_run = "objective")$obj
p <- ob$par
cat("start mu", p[1], " rows with y == start mu:", sum(d$y == p[1]), " n", nrow(d), "\n")
p[1] <- p[1] + 1e-3
cat("AD", signif(as.numeric(ob$gr(p)), 7), "| FD", signif(numDeriv::grad(ob$fn, p), 7), "\n")
