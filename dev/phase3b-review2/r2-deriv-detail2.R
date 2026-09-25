# Reviewer 2: the rows r2-deriv.R flags after its noise floor. Tape
# gradient and Hessian beside Richardson differences at three steps.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(RTMB); library(frmtmb.eam)})
ns <- asNamespace("frmtmb.eam")
options(digits = 10, width = 160)
outs <- list(
  lF = function(p) ns$ddm_rt_lcdf2(p[4], p[1], p[2], p[3])[["lF"]],
  lFl = function(p) ns$ddm_rt_lcdf_b(p[4], p[1], p[2], p[3], 0),
  lFu = function(p) ns$ddm_rt_lcdf_b(p[4], p[1], p[2], p[3], 1),
  lint = function(p) ns$ddm_rt_linterval_b(p[4], p[4] * 1.3, p[1], p[2],
                                           p[3], 0))
rich <- function(f, x, i, h) {
  e <- replace(numeric(length(x)), i, 1)
  D <- function(h) (f(x + h * e) - f(x - h * e)) / (2 * h)
  (4 * D(h / 2) - D(h)) / 3
}
show <- function(nm, p) {
  f <- outs[[nm]]
  tp <- MakeTape(function(q) f(q), p)
  g <- as.numeric(tp$jacobian(p))
  gf <- tp$jacfun(); H <- gf$jacobian(p)
  gg <- function(q) as.numeric(gf(q))
  cat("\n==", nm, paste(signif(p, 6), collapse = " "), " value", f(p), "\n")
  G <- rbind(tape = g)
  for (r in c(1e-4, 1e-5, 1e-6)) {
    G <- rbind(G, vapply(1:4, function(i) rich(f, p, i, r * max(abs(p[i]), 1e-12)), 0))
  }
  rownames(G) <- c("tape", "fd1e-4", "fd1e-5", "fd1e-6"); print(G)
  cat("Hessian: tape, then FD of tape gradient at 1e-5\n")
  print(H)
  print(sapply(1:4, function(i) rich(gg, p, i, 1e-5 * max(abs(p[i]), 1e-12))))
}
show("lF", c(0.7, 1.4, 0.5, 2e-10 * 1.96))
show("lFl", c(0.7, 1.4, 0.5, 2e-10 * 1.96))
show("lint", c(-5.847174087, 2.647920, 0.5609502, 2.133821))
show("lFu", c(0.05 / 1.4, 1.4, 0.3, 2 * 1.96))
show("lFu", c(0.07 / 1.4, 1.4, 0.3, 2 * 1.96))
show("lF", c(1e-9 / 1.4, 1.4, 0.3, 2 * 1.96))
