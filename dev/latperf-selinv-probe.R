# Lane latperf: is a selected inverse of a chain-plus-groups precision
# cheap enough at 32000 rows? TMB's C routine (solveSubset) against an
# R Takahashi recursion. Seed 1.
.libPaths(c("C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(Matrix)
set.seed(1)
mk <- function(n, G) {
  # ar1 chain precision plus G group levels coupled to every row of their
  # group, and a diagonal likelihood weight
  phi <- 0.7
  d <- c(1, rep(1 + phi^2, n - 2), 1)
  Qc <- bandSparse(n, k = c(0, 1), diagonals = list(d + runif(n),
                                                     rep(-phi, n - 1)),
                   symmetric = TRUE)
  g <- sample.int(G, n, TRUE)
  w <- runif(n)
  Z <- sparseMatrix(i = seq_len(n), j = g, x = 1, dims = c(n, G))
  # rows share a latent and a level, so W couples them
  top <- Qc + Diagonal(n, w)
  off <- t(Z) %*% Diagonal(n, w)
  bot <- Diagonal(G, 1) + crossprod(Z * sqrt(w))
  Q <- rbind(cbind(top, t(off)), cbind(off, bot))
  forceSymmetric(as(Q, "CsparseMatrix"))
}
solveSubset <- get("solveSubset", envir = asNamespace("TMB"))
for (n in c(2000, 8000, 32000)) {
  Q <- mk(n, 20)
  t0 <- proc.time()[["elapsed"]]
  L <- Cholesky(Q, super = TRUE, perm = TRUE)
  S <- solveSubset(L = L)
  t1 <- proc.time()[["elapsed"]]
  cat(sprintf("n %d  solveSubset %.2f s  nnz %d\n", n, t1 - t0, length(S@x)))
  if (n <= 2000) {
    D <- as.matrix(solve(Q))
    k <- which(as.matrix(S) != 0, arr.ind = TRUE)
    cat(sprintf("  max rel err on pattern %.2e\n",
                max(abs(S[k] - D[k]) / abs(D[k]))))
  }
}
