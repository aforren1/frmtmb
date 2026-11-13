# Lane arspat, fix round (review MINOR-2): rho's interval from the real
# eigenvalues alone (spatialreg's rule), and log|I - rho W| against
# determinant() inside it and close to both ends, on a W whose complex
# real part (-0.5) is more extreme than its real eigenvalues (+/- 0.2),
# and on the lane's 3-NN test W. Output: dev/arspat-out/sar-ends.txt
.libPaths(c("C:/Users/adf44/source/r/wt-arspat-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.spatial))
blk <- matrix(0, 5, 5)
blk[cbind(1:3, c(2, 3, 1))] <- 1
blk[4, 5] <- blk[5, 4] <- 0.2
Wb <- kronecker(diag(4), blk)
set.seed(21)
n <- 40
xy <- cbind(runif(n), runif(n))
D <- as.matrix(dist(xy)); diag(D) <- Inf
K <- t(apply(D, 1, function(r) as.numeric(rank(r) <= 3)))
Wk <- K / rowSums(K)
for (nm in c("block W", "3-NN W")) {
  W <- if (nm == "block W") Wb else Wk
  nn <- nrow(W)
  ev <- eigen(W, only.values = TRUE)$values
  real <- abs(Im(ev)) <= 1e-8 * max(Mod(ev))
  aux <- frmtmb.spatial:::sar_build(list(M = quote(W), type = "lag"),
                                    seq_len(nn), nn, function(e, w) W)
  cat(sprintf("%s: interval (%.6f, %.6f); 1/range(real eig) (%.6f, %.6f); 1/range(Re all) (%.6f, %.6f)\n",
              nm, aux$lo, aux$hi, 1 / min(Re(ev[real])), 1 / max(Re(ev[real])),
              1 / min(Re(ev)), 1 / max(Re(ev))))
  pts <- c(aux$lo * (1 - 1e-6), aux$lo / 2, 0, aux$hi / 2,
           aux$hi * (1 - 1e-6))
  for (rho in pts) {
    dt <- determinant(diag(nn) - rho * W)
    ld <- as.numeric(frmtmb.spatial:::sar_logdet(rho, aux))
    cat(sprintf("  rho %+.6f: sign %+d, log|det| %.6f, |diff| %.2e\n", rho,
                as.integer(dt$sign), as.numeric(dt$modulus),
                abs(ld - as.numeric(dt$modulus))))
  }
}
