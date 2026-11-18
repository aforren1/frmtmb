# Short check of lane sparsese: jcs_fd_root() solves Ld X = V_db with
# qr.solve(), whose default tol = 1e-7 declares a factor rank-deficient
# when the drawn parameters' standard deviations differ by about 1e7 and
# they correlate. Does it stop where jcs_psd_factor() succeeded?
#   Rscript sparsese-short-qrsolve.R
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
pf <- frmtmb:::jcs_psd_factor
for (ratio in 10^c(3, 5, 6, 7, 8, 9)) {
  for (rho in c(0.5, 0.99)) {
    s <- c(1, 2, ratio)
    Cr <- matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), 3)
    Vdd <- Cr * outer(s, s)
    Vdb <- matrix(c(0.3, 0.2, 0.1) * s, 3, 1)  # a group effect, sd 1
    Ld <- pf(Vdd)
    r <- tryCatch({
      Kb <- t(qr.solve(Ld, Vdb))
      sprintf("ok, |Kb Ld' - V_bd| rel %.1e", max(abs(Kb %*% t(Ld) - t(Vdb))) /
                max(abs(Vdb)))
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("sd ratio %.0e rho %.2f: %s\n", ratio, rho, r))
  }
}
# the same solve without qr's rank tolerance: the factor is full column
# rank by construction (a Cholesky factor, or eigen-directions kept above
# 1e-12)
for (ratio in 10^c(8, 9, 12)) {
  s <- c(1, 2, ratio)
  Cr <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), 3)
  Vdd <- Cr * outer(s, s)
  Vdb <- matrix(c(0.3, 0.2, 0.1) * s, 3, 1)
  Ld <- pf(Vdd)
  Kb <- t(qr.solve(Ld, Vdb, tol = 0))
  Vc <- 1 - t(Vdb / s) %*% solve(Cr, Vdb / s)
  cat(sprintf(paste0("tol = 0, sd ratio %.0e: |Kb Ld' - V_bd| rel %.1e;",
                     " 1 - Kb Kb' %.6f vs conditional variance %.6f\n"),
              ratio, max(abs(Kb %*% t(Ld) - t(Vdb))) / max(abs(Vdb)),
              1 - sum(Kb^2), Vc))
}
