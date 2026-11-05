# Does the correlation plus a relative 1e-12 ever fail to factor in double?
# Seed 20261103. Designs: 1-D uniform and grid, 2-D, duplicates at gr = FALSE,
# clustered positions, n to 500. Scans the length scale on a fine grid, also at
# smaller nuggets, so the margin is a number. Uses R's LAPACK.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
set.seed(20261103)
des <- list()
for (n in c(60L, 200L, 500L)) {
  des[[paste0("unif1d_", n)]] <- matrix(runif(n, 0, 10))
  des[[paste0("grid1d_", n)]] <- matrix(seq(0, 10, length.out = n))
  des[[paste0("cluster_", n)]] <- matrix(c(runif(n - 5, 0, 1e-3), runif(5, 0, 10)))
  des[[paste0("unif2d_", n)]] <- cbind(runif(n, 0, 10), runif(n, 0, 10))
  u <- runif(ceiling(n / 10), 0, 10)
  des[[paste0("dup10_", n)]] <- matrix(rep(u, length.out = n))
  des[[paste0("dupall_", n)]] <- matrix(rep(1.5, n) + c(rep(0, n - 2), 1, 2))
}
grid <- exp(seq(log(1e-3), log(1e7), length.out = 240))
# the order frm() factors in; dev/ncgp-perm-patch.R's and the frame's
far_first <- function(d2) {
  m <- nrow(d2); p <- integer(m); p[1L] <- which.min(rowSums(d2))
  gap <- d2[, p[1L]]; gap[p[1L]] <- -1
  for (k in seq_len(m - 1L) + 1L) {
    p[k] <- which.max(gap); gap <- pmin(gap, d2[, p[k]]); gap[p[seq_len(k)]] <- -1
  }
  p
}
for (dn in names(des)) {
  X <- des[[dn]]
  D2 <- as.matrix(dist(X))^2
  po <- far_first(D2)
  D2 <- D2[po, po]
  for (nug in c(1e-12, 1e-13, 1e-14, 1e-15, 1e-16)) {
    fails <- 0L
    first <- NA
    for (r in grid) {
      K <- exp(-D2 / (2 * r^2))
      diag(K) <- diag(K) + nug
      ok <- !inherits(try(chol(K), silent = TRUE), "try-error")
      if (!ok) {
        fails <- fails + 1L
        if (is.na(first)) first <- r
      }
    }
    cat(sprintf("NCGPS %-12s nugget %.0e: fails %3d of %d, first rho %s\n", dn, nug,
                fails, length(grid), format(first, digits = 3)))
    if (fails == 0L && nug < 1e-14) break
  }
}
