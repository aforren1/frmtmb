# Punch 1, B2: does the value and the gradient join across every switch?
#   Rscript dev/fams2-p1-joins.R <impl>   ("cand" or "pkg")
# At each edge (x = m, x = m +- sd, |x - m| = 5 sd and 6 sd, a b / (a + b) =
# 150 and 450), the function is read at the edge +- a relative 1e-9 and
# the jump is compared with what the gradient itself predicts over that
# step. Then the reviewer's point a = 3e4, b = 7e4 across x = m.
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
f <- if (args[1] == "cand") {
  source("C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-p1-cand.R")
  cand
} else frmtmb:::log_ibeta_half
F <- MakeTape(function(p) f(p[1], p[2], p[3]), c(0.2, 3, 5))
J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
jump <- function(p0, dir, lab) {
  # a relative step in the coordinate that moves
  e <- 1e-9 * abs(sum(p0 * dir)) * dir
  lo <- p0 - e; hi <- p0 + e
  dv <- F(hi) - F(lo)
  pred <- sum((J(lo) + J(hi)) / 2 * (hi - lo))
  gj <- max(abs(J(hi) - J(lo))) / max(1, abs(J(p0)))
  cat(sprintf("  %-34s value jump - predicted %9.2e  gradient jump (rel) %9.2e  finite H,T3 at edge %s\n",
              lab, dv - pred, gj,
              all(is.finite(H(p0))) && all(is.finite(T3(p0)))))
}
for (ab in list(c(0.4, 2), c(5, 40), c(300, 700), c(3e4, 7e4), c(2e5, 8e5))) {
  a <- ab[1]; b <- ab[2]; s <- a + b
  m <- (a + 1) / (s + 2); sd <- sqrt(a * b / (s * s * (s + 1)))
  cat(sprintf("a = %g, b = %g (m = %.6g, sd = %.3g)\n", a, b, m, sd))

  for (kk in c(0, -1, 1, 5, -5, 6, -6)) {
    x <- m + kk * sd
    if (x <= 0 || x >= 0.5) next
    jump(c(x, a, b), c(1, 0, 0), sprintf("x = m %+g sd", kk))
  }
}
cat("across min(a, b) = 150 and 450, x at m + 0.5 sd and m + 3.5 sd:\n")
for (mn in c(150, 450)) for (kk in c(0.5, 5.5)) {
  a <- mn / 0.8; b <- 4 * a; s <- a + b
  m <- (a + 1) / (s + 2); sd <- sqrt(a * b / (s * s * (s + 1)))
  jump(c(m + kk * sd, a, b), c(0, 1, 0), sprintf("a = %g, k = %g", mn, kk))
}
cat("the reviewer's point, a = 3e4, b = 7e4, either side of x = m:\n")
a <- 3e4; b <- 7e4; m <- (a + 1) / (a + b + 2)
for (e in c(-1e-12, 0, 1e-12)) {
  x <- m + e
  g <- J(c(x, a, b))
  ref <- exp(dbeta(x, a, b, log = TRUE) - pbeta(x, a, b, log.p = TRUE))
  cat(sprintf("  x = m %+.0e: d/dx %.6f  reference %.6f  value %.15f  pbeta %.15f\n",
              e, g[1], ref, F(c(x, a, b)), pbeta(x, a, b, log.p = TRUE)))
}
