# Reviewer 2: detail of the worst rows of r2-deriv.R. Prints the tape
# gradient beside Richardson differences at two step sizes, so a bad
# finite difference can be told from a bad tape.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(RTMB); library(frmtmb.eam)})
ns <- asNamespace("frmtmb.eam")
options(digits = 12, width = 160)
f_lF <- function(p) ns$ddm_rt_lcdf2(p[4], p[1], p[2], p[3])[["lF"]]
f_lFl <- function(p) ns$ddm_rt_lcdf_b(p[4], p[1], p[2], p[3], 0)
rich <- function(f, x, i, h) {
  e <- replace(numeric(length(x)), i, 1)
  D <- function(h) (f(x + h * e) - f(x - h * e)) / (2 * h)
  (4 * D(h / 2) - D(h)) / 3
}
show <- function(f, p, lab) {
  tp <- MakeTape(function(q) f(q), p)
  g <- as.numeric(tp$jacobian(p))
  fd5 <- vapply(1:4, function(i) rich(f, p, i, 1e-5 * max(abs(p[i]), 1e-3)), 0)
  fd3 <- vapply(1:4, function(i) rich(f, p, i, 1e-3 * max(abs(p[i]), 1e-3)), 0)
  cat("\n==", lab, " value", f(p), "\n")
  print(rbind(tape = g, fd_1e5 = fd5, fd_1e3 = fd3))
}
show(f_lF, c(0, 1.4, 0.4, 9.8), "lF v=0 u=5")
show(f_lF, c(0.7, 1.4, 0.5, 0.05 * exp(4.8 - 0.02) * 1.96), "lF clamp S +")
show(f_lFl, c(0.7, 1.4, 0.5, 1.001e-10 * 1.96), "lFl hold 1.001e-10")
show(f_lFl, c(0.7, 1.4, 0.5, 2e-10 * 1.96), "lFl hold 2e-10")
show(f_lFl, c(0.7, 1.4, 0.5, 1e-6 * 1.96), "lFl u=1e-6")
show(f_lFl, c(0.7, 1.4, 0.5, 1e-4 * 1.96), "lFl u=1e-4")
# analytic leading term at small u: log F_l ~ -w^2 a^2 / (2 t) + ...
p <- c(0.7, 1.4, 0.5, 1.001e-10 * 1.96)
cat("\nleading-term d/dt of -z^2/(2t):", (0.7^2 * 0 + (0.5 * 1.4)^2) / (2 * p[4]^2), "\n")
