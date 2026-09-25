# Second derivatives: the Laplace approximation needs the Hessian in
# the random effects, so a NaN there breaks a fit the first-derivative
# hunt (phase3b-nan-hunt4.R) passes.
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
  logFb = function(p) sum(log(ns$ddm_floor(exp(ns$ddm_rt_lcdf_b(q - p[3], p[1],
                                                        exp(p[2]), 0.5,
                                                        up)), 1e-300))),
  dens = function(p) sum(ns$ddm_lpdf_both(q - p[3], p[1], exp(p[2]), 0.5, up)))
grid <- expand.grid(v = c(-20, -5, -1, 0, 1, 5, 20),
                    la = c(-4, -2, -0.5, 0, 0.5, 1.5, 3),
                    t0 = c(0, 0.1, 0.25, 0.4, 0.4499, 0.45))
for (nm in names(pieces)) {
  tp <- RTMB::MakeTape(pieces[[nm]], c(0, 0, 0))
  hf <- tp$jacfun()
  bad <- 0; first <- NULL
  for (i in seq_len(nrow(grid))) {
    p <- as.numeric(grid[i, ])
    h <- hf$jacobian(p)
    if (!all(is.finite(h))) {
      bad <- bad + 1
      if (is.null(first)) first <- p
    }
  }
  cat(nm, ": non-finite Hessian at", bad, "of", nrow(grid))
  if (!is.null(first)) cat("; first:", format(first, digits = 4))
  cat("\n")
}
# where in the survival
sub <- list(
  small = function(p) sum(ns$ddm_lsurv_small(ns$ddm_floor(q - p[3], 1e-10 * exp(2 * p[2])), p[1], exp(p[2]), 0.5)),
  large = function(p) sum(ns$ddm_lsurv_large(ns$ddm_floor(q - p[3], 1e-10 * exp(2 * p[2])), p[1], exp(p[2]), 0.5)),
  lpm = function(p) sum(ns$ddm_lpmass(p[1] * q - 3, p[1] * q + p[2])))
for (nm in names(sub)) {
  tp <- RTMB::MakeTape(sub[[nm]], c(0, 0, 0)); hf <- tp$jacfun()
  bad <- sum(apply(grid, 1, function(p) !all(is.finite(hf$jacobian(as.numeric(p))))))
  cat(nm, ": non-finite Hessian at", bad, "of", nrow(grid), "\n")
}
