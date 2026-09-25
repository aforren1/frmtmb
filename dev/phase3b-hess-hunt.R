# Non-finite second derivatives of the distribution-function pieces
# over ordinary values: drift -4..4, boundary 0.5..4, decision times
# 0.02..3.
a0 <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a0[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb.eam")
tt <- c(0.02, 0.1, 0.225, 0.5, 1, 3)
pieces <- list(
  llower_small = function(p, t) ns$ddm_llower_small(t, p[1], exp(p[2]), 0.5),
  llower_large = function(p, t) ns$ddm_llower_large(t, p[1], exp(p[2]), 0.5),
  ltail = function(p, t) ns$ddm_ltail_lower_large(t, p[1], exp(p[2]), 0.5),
  lprob = function(p, t) ns$ddm_llower_prob_s(p[1], exp(p[2]), 0.5) + 0 * t,
  lcdf_b = function(p, t) ns$ddm_rt_lcdf_b(t, p[1], exp(p[2]), 0.5, 1),
  linterval = function(p, t) ns$ddm_rt_linterval_b(t, t + 0.1, p[1], exp(p[2]), 0.5, 1),
  lF = function(p, t) ns$ddm_rt_lcdf2(t, p[1], exp(p[2]), 0.5)$lF)
grid <- expand.grid(v = c(-4, -1.5, -0.3, 0, 0.02, 0.3, 1.5, 4),
                    la = log(c(0.5, 1, 1.5, 2.5, 4)))
for (nm in names(pieces)) {
  bad <- NULL
  for (t in tt) {
    tp <- RTMB::MakeTape(function(p) pieces[[nm]](p, t), c(0, 0))
    hf <- tp$jacfun()
    for (i in seq_len(nrow(grid))) {
      p <- as.numeric(grid[i, ])
      if (!all(is.finite(hf$jacobian(p)))) bad <- rbind(bad, c(t, p))
    }
  }
  cat(nm, ": non-finite Hessian at", NROW(bad), "of", length(tt) * nrow(grid))
  if (!is.null(bad)) cat("; e.g. t, v, log a =", format(bad[1, ], digits = 3))
  cat("\n")
}
