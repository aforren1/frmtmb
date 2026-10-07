# Lane optima, item 2: how far into the tail is RTMB's pnorm(log.p =
# TRUE) accurate, in value and derivative, against Rmpfr at 256 bits?
# x on a log grid; the lower tail log Phi(-x) and its derivative
# -phi(x) / Phi(-x).
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
x <- 10^seq(0.5, 8, by = 0.05)
tp <- MakeTape(function(e) RTMB::pnorm(-e, log.p = TRUE), x)
v <- tp(x)
g <- diag(tp$jacobian(x))
X <- mpfr(x, 256)
P <- Rmpfr::pnorm(-X)
rv <- as.numeric(log(P))
rg <- as.numeric(-Rmpfr::dnorm(X) / P)
ev <- abs(v / rv - 1)
eg <- abs(g / rg - 1)
for (cut in c(1e2, 1e3, 1e4, 3e4, 1e5, 3e5, 1e6, 1e7, 1e8)) {
  k <- x <= cut
  cat(sprintf("x <= %-6g max rel err value %.3g, derivative %.3g\n", cut,
              max(ev[k]), max(eg[k])))
}
