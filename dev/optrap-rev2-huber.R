# Reviewer of lane optrap, re-check of m6: is huber()'s Laplace objective
# discontinuous on seed 19 (dev/optrap-p1-huber19.R's segment)? Take the
# largest step on a 2001-point grid of the segment and bisect it: a jump
# keeps its size as the interval shrinks; a steep but smooth stretch
# does not.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(19)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
b <- rnorm(20)[d$g]
d$y <- (1 + 0.5 * d$x + 0.5 * b + rt(100, 3)) * 1000
f <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = huber(), data = d))
pp <- f$opt$par
p7 <- readRDS("dev/optrap-log/p1-huber19-r7.rds")
p7[1:2] <- p7[1:2] * 1000
p7[3:4] <- p7[3:4] + log(1000)
fn <- function(s) -f$obj$fn(pp + s * (p7 - pp))
t <- seq(0, 1, length.out = 2001)
v <- vapply(t, fn, 0)
dv <- diff(v)
cat("lane logLik", fn(0), "r7 point", fn(1), "| max on grid", max(v), "\n")
for (k in order(-abs(dv))[1:3]) {
  lo <- t[k]; hi <- t[k + 1]
  cat(sprintf("grid step at t = %.4f: %.4g\n", lo, dv[k]))
  for (it in 1:30) {
    mid <- (lo + hi) / 2
    if (abs(fn(mid) - fn(lo)) >= abs(fn(hi) - fn(mid))) hi <- mid else lo <- mid
    if (it %% 10 == 0) {
      cat(sprintf("  interval %.3g: change across it %.4g\n", hi - lo,
                  fn(hi) - fn(lo)))
    }
  }
}
