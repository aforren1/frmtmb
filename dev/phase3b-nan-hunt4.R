# NaN hunt on the rewritten distribution functions: tape each and scan
# parameter values an optimizer step can reach, including a decision
# time at or below zero.
a0 <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a0[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb.eam")
q <- c(0.45, 0.5, 0.6, 1.0, 2.5)
up <- c(0, 1, 1, 0, 1)
pieces <- list(
  lS = function(p) sum(ns$ddm_rt_lcdf2(q - p[3], p[1], exp(p[2]), 0.5)$lS),
  lF = function(p) sum(ns$ddm_rt_lcdf2(q - p[3], p[1], exp(p[2]), 0.5)$lF),
  lFb = function(p) sum(ns$ddm_rt_lcdf_b(q - p[3], p[1], exp(p[2]), 0.5, up)),
  Fb = function(p) sum(exp(ns$ddm_rt_lcdf_b(q - p[3], p[1], exp(p[2]), 0.5,
                                            up))))
grid <- expand.grid(v = c(-60, -20, -5, -1, 0, 1e-9, 1, 5, 20, 60),
                    la = c(-8, -4, -2, -0.5, 0, 0.5, 1.5, 3, 6),
                    t0 = c(-5, 0, 0.1, 0.25, 0.4499999, 0.45, 0.5, 2))
for (nm in names(pieces)) {
  tp <- RTMB::MakeTape(pieces[[nm]], c(0, 0, 0))
  bad <- 0; first <- NULL
  for (i in seq_len(nrow(grid))) {
    p <- as.numeric(grid[i, ])
    val <- tp(p); g <- tp$jacobian(p)
    if (!all(is.finite(g)) || !is.finite(val)) {
      bad <- bad + 1
      if (is.null(first)) first <- c(p, val, g)
    }
  }
  cat(nm, ": non-finite at", bad, "of", nrow(grid))
  if (!is.null(first)) cat("; first:", format(first, digits = 4))
  cat("\n")
}
