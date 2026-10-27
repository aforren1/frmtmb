# Reviewer of lane loglik, fix round 1: why a direct call of the
# jacobian check fails on a 1600-row sar() that frm() accepted.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
set.seed(101)
for (k in c(6, 40)) {
  nn <- k * k
  id <- matrix(seq_len(nn), k, k)
  Pn <- matrix(0, nn, nn)
  Pn[cbind(as.vector(id[-k, ]), as.vector(id[-1, ]))] <- 1
  Pn[cbind(as.vector(id[, -k]), as.vector(id[, -1]))] <- 1
  Pn <- Pn + t(Pn)
  W16 <- Pn / rowSums(Pn)
  d16 <- data.frame(x = rnorm(nn))
  d16$y <- 1 + rnorm(nn)
  f16 <- frm(y ~ x + sar(W16), data = d16, data2 = list(W16 = W16),
             dry_run = "objective")
  ac <- f16$frame$autocor$y
  e <- frmtmb:::ac_block_entry(ac)
  st <- e$start(ac$aux)
  i <- seq_len(nn)
  dirs <- cbind(sin(i), cos(1.7 * i))
  m <- as.numeric(e$unwhiten(numeric(nn), numeric(nn), 1, st, ac$aux))
  Y <- e$unwhiten(dirs, numeric(nn), 1, st, ac$aux)
  cat("k", k, "class of unwhiten(matrix):", class(Y), "dim", dim(Y), "\n")
  Y <- as.matrix(Y)
  U <- vapply(1:2, function(j) as.numeric(e$whiten(Y[, j], numeric(nn), 1, st,
                                                    ac$aux)$u), numeric(nn))
  K <- e$jacobian(1, st, ac$aux)
  KD <- as.matrix(K %*% (Y - m))
  cat("  max|KD - U| =", max(abs(KD - U)), " max|U| =", max(abs(U)),
      " max|U - dirs| =", max(abs(U - dirs)), "\n")
  r <- tryCatch({
    frmtmb:::autocor_seam_check_jacobian(e, ac$aux, st, nn,
                                         function(w) stop(w))
    "passes"
  }, error = function(err) substr(conditionMessage(err), 1, 200))
  cat("  direct check:", r, "\n")
}
