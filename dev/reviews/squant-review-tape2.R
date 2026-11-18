# Review of lane squant: locate the zero tape gradient of lF in
# quantile_elf()'s cdf_log (squant-review-tape.txt (i)).
# Run (worktree root):
#   Rscript dev/reviews/squant-review-tape2.R > dev/reviews/squant-review-tape2.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
a <- 0.125
b <- 0.125
jac <- function(f, x0) as.numeric(RTMB::MakeTape(f, x0)$jacobian(x0))
cat("t, d lF/dt (cdf_log via t), d log_ibeta_two$lF/dt, d/dt of x, d log_ibeta_half/dx * x(1-x), true a\n")
for (t in c(-50, -100, -150, -200, -250, -300, -350, 50, 300)) {
  g1 <- jac(function(tt) {
    ns$elf_cdf_log(tt, 0, 1, 0.5, 0.5, 0.25)[[if (t < 0) "lF" else "lS"]]
  }, t * 0.25)
  g2 <- jac(function(tt) {
    x <- 1 / (1 + exp(-tt))
    xc <- 1 / (1 + exp(tt))
    ns$log_ibeta_two(x, xc, a, b)[[if (t < 0) "lF" else "lS"]]
  }, t)
  g3 <- jac(function(tt) 1 / (1 + exp(-tt)), t)
  x0 <- 1 / (1 + exp(-t))
  g4 <- jac(function(x) ns$log_ibeta_half(x, a, b), x0) * x0 * (1 - x0)
  g5 <- jac(function(tt) {
    x <- 1 / (1 + exp(-tt))
    xc <- 1 / (1 + exp(tt))
    ns$log_ibeta_two(x, xc, a, b)[["lS"]]
  }, t)
  cat(sprintf("%5g: %.6g (per unit t: %.6g) | %.6g | %.6g | %.6g | lS %.6g\n",
              t, g1, g1 * 0.25, g2, g3, g4, g5))
}
cat("DONE\n")
