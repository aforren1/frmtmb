# Short check of lane sparsese, fix round 1: jcs_psd_factor() on PSD inputs.
#   Rscript sparsese-short-psd.R
# For each input V: max |L L' - V| / max |V|, whether any diagonal entry of
# L L' exceeds V's (invented variance), the rank kept, and which branch.
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
pf <- frmtmb:::jcs_psd_factor
set.seed(11)
rnd <- function(n, k, sc = rep(1, n)) {
  B <- matrix(stats::rnorm(n * k), n, k) * sc
  tcrossprod(B)
}
chk <- function(nm, V) {
  L <- pf(V)
  R <- tcrossprod(L)
  br <- if (ncol(L) == sum(diag(V) > 0)) "chol(or eig full)" else "eigen"
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
  dd <- diag(R) - diag(V)
  zero <- which(diag(V) == 0)
  cat(sprintf(paste0("%-34s n=%3d rank(V, 1e-12)=%3d cols=%3d  ",
                     "max|LL'-V|/max|V| %.2e  max(diag LL' - diag V)/max",
                     "|V| %+.2e  zero rows of L at zero-var: %s\n"),
              nm, nrow(V), sum(ev > max(ev) * 1e-12), ncol(L),
              max(abs(R - V)) / max(abs(V)), max(dd) / max(abs(V)),
              if (length(zero)) all(L[zero, ] == 0) else "n/a"))
  # entrywise in correlation scale, where a block of huge and tiny
  # standard deviations is judged on each pair, not on the largest
  s <- sqrt(pmax(diag(V), 0))
  ok <- s > 0
  cat(sprintf("%-34s scaled max|LL'-V|/(s_i s_j) %.2e\n", "",
              max(abs(R - V)[ok, ok] / outer(s[ok], s[ok]))))
}
# 1. full rank, widely scaled (natural-coordinate scales 1e-4 to 1e4)
V1 <- rnd(8, 12, 10^seq(-4, 4, length.out = 8))
chk("full rank, scales 1e-4..1e4", V1)
# 2. zero rows and columns (a lost parameter of a repaired covariance)
V2 <- matrix(0, 10, 10)
V2[-c(3, 7), -c(3, 7)] <- rnd(8, 12)
chk("zero rows 3 and 7", V2)
# 3. rank-deficient block, exact (rank 4 of 9)
V3 <- rnd(9, 4)
chk("rank 4 of 9", V3)
# 4. rank-deficient plus zero rows
V4 <- matrix(0, 11, 11)
V4[-c(1, 11), -c(1, 11)] <- rnd(9, 5)
chk("rank 5 of 9, zero rows 1 and 11", V4)
# 5. nearly singular: one eigenvalue 1e-10, 1e-14, 1e-18 of the largest
for (tiny in c(1e-10, 1e-14, 1e-18)) {
  Q <- qr.Q(qr(matrix(stats::rnorm(64), 8, 8)))
  V5 <- Q %*% diag(c(1, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01, tiny)) %*% t(Q)
  V5 <- (V5 + t(V5)) / 2
  chk(sprintf("near singular, min eig %.0e", tiny), V5)
}
# 6. near-singular from perfectly correlated pair, scaled apart
V6 <- rnd(6, 10)
V6 <- rbind(cbind(V6, V6[, 1] * 1e3), c(V6[1, ] * 1e3, V6[1, 1] * 1e6))
chk("col 7 = 1e3 * col 1 (rank 6 of 7)", V6)
# 7. a tiny negative eigenvalue from rounding (indefinite by 1e-16)
V7 <- rnd(6, 3)
V7 <- V7 - 1e-16 * max(abs(V7)) * diag(6)
chk("rank 3 of 6, shifted -1e-16", V7)
# 8. all zero
L8 <- pf(matrix(0, 3, 3))
cat("all zero: dim(L)", dim(L8), "\n")
# 9. a 1 x 1
cat("1x1 (4): L =", pf(matrix(4, 1, 1)), "\n")
