# Second derivatives of the distribution-function pieces against
# central differences of the tape gradient, at ordinary values.
a0 <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a0[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb.eam")
t <- c(0.05, 0.225, 0.3, 0.8, 2)
pieces <- list(
  llower_small = function(p) sum(ns$ddm_llower_small(t, p[1], exp(p[2]), 0.5)),
  llower_large = function(p) sum(ns$ddm_llower_large(t, p[1], exp(p[2]), 0.5)),
  lcdf_b = function(p) sum(ns$ddm_rt_lcdf_b(t, p[1], exp(p[2]), 0.5, c(0, 1, 0, 1, 0))),
  linterval = function(p) sum(ns$ddm_rt_linterval_b(t, t + 0.1, p[1], exp(p[2]), 0.5, c(0, 1, 0, 1, 0))),
  lS = function(p) sum(ns$ddm_rt_lcdf2(t, p[1], exp(p[2]), 0.5)$lS),
  lF = function(p) sum(ns$ddm_rt_lcdf2(t, p[1], exp(p[2]), 0.5)$lF),
  lpnorm = function(p) sum(ns$ddm_lpnorm(p[1] * c(-60, -10, -1, 2, 40) + p[2])),
  lpmass = function(p) sum(ns$ddm_lpmass(p[1] * c(-3, 1, 5) + p[2], p[1] * c(-1, 2, 9) + p[2] + 1)),
  dens = function(p) sum(ns$ddm_lpdf_both(t, p[1], exp(p[2]), 0.5, c(0, 1, 0, 1, 0))))
for (p in list(c(0, log(1.5)), c(1.3, log(1.5)), c(-2, log(0.8)))) {
  for (nm in names(pieces)) {
    tp <- RTMB::MakeTape(pieces[[nm]], p)
    H <- tp$jacfun()$jacobian(p)
    fd <- sapply(1:2, function(k) {
      h <- 1e-5; e <- replace(numeric(2), k, h)
      (tp$jacobian(p + e) - tp$jacobian(p - e)) / (2 * h)
    })
    cat(sprintf("p=(%g,%.3g) %-13s max|H - fd| / max(1,|fd|) = %.3g   H[1,1] %.4g fd %.4g\n",
                p[1], p[2], nm, max(abs(H - fd) / pmax(1, abs(fd))), H[1, 1], fd[1, 1]))
  }
}
