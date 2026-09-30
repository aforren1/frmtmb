# Punch 1, B2: sweep a log incomplete beta implementation.
#   Rscript dev/fams2-p1-sweep.R <impl> [points]
# <impl> is "cand" (dev/fams2-p1-cand.R) or "pkg" (the installed
# frmtmb:::log_ibeta_half()). Seed 20260930. Shapes log-uniform on
# [1e-3, 1e7] each; half the x drawn within 5 sd of m, a tenth exactly at
# m, the rest uniform on (0, 1/2). Value error against stats::pbeta(log.p
# = TRUE), floored at one (the reviewer found stats::pbeta() and Rmpfr
# to agree to 4e-14 here); a subsample for the gradient.
args <- commandArgs(TRUE)
impl <- args[1]
M <- if (length(args) >= 2) as.integer(args[2]) else 4000L
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
f <- if (impl == "cand") {
  source("C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-p1-cand.R")
  cand
} else {
  frmtmb:::log_ibeta_half
}
set.seed(20260930)
a <- exp(runif(M, log(1e-3), log(1e7)))
b <- exp(runif(M, log(1e-3), log(1e7)))
s <- a + b
m <- (a + 1) / (s + 2)
sd <- sqrt(a * b / (s * s * (s + 1)))
u <- runif(M)
x <- ifelse(u < 0.1, m, ifelse(u < 0.6, m + runif(M, -5, 5) * sd,
                                runif(M, 0, 0.5)))
ok <- x > 0 & x < 0.5
x <- x[ok]; a <- a[ok]; b <- b[ok]; m <- m[ok]; sd <- sd[ok]
ref <- pbeta(x, a, b, log.p = TRUE)
keep <- is.finite(ref) & ref > -600
x <- x[keep]; a <- a[keep]; b <- b[keep]; ref <- ref[keep]
m <- m[keep]; sd <- sd[keep]
n <- length(x)
F <- MakeTape(function(p) f(p[1], exp(p[2]), exp(p[3])), c(0.2, 0, 0))
J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
val <- numeric(n); fin <- matrix(NA, n, 3)
for (i in seq_len(n)) {
  p <- c(x[i], log(a[i]), log(b[i]))
  val[i] <- F(p)
  fin[i, ] <- c(all(is.finite(J(p))), all(is.finite(H(p))),
                all(is.finite(T3(p))))
}
err <- abs(val - ref) / pmax(1, abs(ref))
err[!is.finite(err)] <- Inf
k <- abs(x - m) / sd
mn <- pmin(a, b)
cat(impl, ": points", n, "(of", M, "drawn; log P > -600 and 0 < x < 1/2)\n")
cat("non-finite: value", sum(!is.finite(val)), " gradient", sum(!fin[, 1]),
    " hessian", sum(!fin[, 2]), " third", sum(!fin[, 3]), "\n")
cat(sprintf("value error: max %.2e  q99 %.2e  q999 %.2e\n", max(err),
            quantile(err, 0.99), quantile(err, 0.999)))
options(width = 200)
cat("log10 worst value error, |x - m| / sd (rows) by max(a, b) (cols):\n")
print(round(log10(pmax(tapply(err, list(
  cut(k, c(-1, 1e-9, 0.1, 1, 3, 5, Inf)),
  cut(pmax(a, b), c(0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e8))), max), 1e-17)), 1))
cat("at x == m exactly:", sum(k < 1e-12), "points, worst",
    format(max(c(0, err[k < 1e-12])), digits = 3), "\n")
# the gradient in (x, log a, log b): the x derivative exactly,
# dbeta / I, and the shape derivatives by a five-point stencil of
# stats::pbeta(log.p = TRUE), h = min(1e-3, 0.02 / sqrt(a + b)) so the
# stencil moves m by a small part of a standard deviation
set.seed(7)
idx <- sample(n, min(400L, n))
lref <- function(xx, la, lb) pbeta(xx, exp(la), exp(lb), log.p = TRUE)
st <- function(fn, h) (-fn(2 * h) + 8 * fn(h) - 8 * fn(-h) + fn(-2 * h)) / (12 * h)
ge <- numeric(length(idx))
for (j in seq_along(idx)) {
  i <- idx[j]
  p <- c(x[i], log(a[i]), log(b[i]))
  gx <- exp(dbeta(x[i], a[i], b[i], log = TRUE) - ref[i])
  hh <- min(1e-3, 0.02 / sqrt(a[i] + b[i]))
  ga <- st(function(h) lref(p[1], p[2] + h, p[3]), hh)
  gb <- st(function(h) lref(p[1], p[2], p[3] + h), hh)
  gr <- c(gx, ga, gb)
  g <- J(p)
  ge[j] <- max(abs(g - gr) / pmax(1, abs(gr)))
}
cat(sprintf("gradient against the reference (%d points): max %.2e  q99 %.2e  median %.2e\n",
            length(idx), max(ge), quantile(ge, 0.99), median(ge)))
print(summary(data.frame(a = a[idx], b = b[idx], k = (x[idx] - m[idx]) / sd[idx], ge)[order(-ge)[1:5], ]))
# the worst gradient points against RTMB::pbeta()'s own AD gradient,
# which is exact wherever it is finite
Fp <- MakeTape(function(p) log(RTMB::pbeta(p[1], exp(p[2]), exp(p[3]))),
               c(0.2, 0, 0))
Jp <- Fp$jacfun()
worst <- idx[order(-ge)[1:10]]
gp <- t(vapply(worst, function(i) {
  p <- c(x[i], log(a[i]), log(b[i]))
  g <- J(p); r <- Jp(p)
  c(ours_vs_rtmb = max(abs(g - r) / pmax(1, abs(r))),
    stencil_vs_rtmb = ge[match(i, idx)])
}, numeric(2)))
cat("the ten worst: ours against RTMB AD, max", format(max(gp[, 1]), digits = 3),
    "; the stencil's own disagreement there, max", format(max(gp[, 2]), digits = 3), "\n")
