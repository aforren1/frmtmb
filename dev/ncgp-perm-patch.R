# The farthest-first arm: frm()'s factor over gpgap's order of the positions, patched
# into this process's frmtmb namespace. The kernel and the model are the same; the
# rounding and the basis of z move. Source it after library(frmtmb).
#   maximin_D2() as dev/review-gpgap-perm.R in the gpgap tree.
ncgp_maximin <- function(D2) {
  m <- nrow(D2)
  p <- which.min(rowSums(D2))
  while (length(p) < m) {
    d <- apply(D2[, p, drop = FALSE], 1L, min)
    p <- c(p, which.max(d))
  }
  p
}
ncgp_order_cache <- new.env()
ncgp_perm <- function(blk) {
  D2 <- Reduce(`+`, blk[["aux_D2"]])
  key <- paste(nrow(D2), sum(D2))
  if (is.null(ncgp_order_cache[[key]])) ncgp_order_cache[[key]] <- ncgp_maximin(D2)
  ncgp_order_cache[[key]]
}
local({
  ns <- asNamespace("frmtmb")
  gp_corr <- get("gp_corr", ns)
  gp_chol_num <- get("gp_chol_num", ns)
  sd_floored <- get("sd_floored", ns)
  gp_nugget <- get("gp_nugget", ns)
  perm_factor <- function(theta, blk) {
    K <- gp_corr(theta, blk)
    p <- ncgp_perm(blk)
    Kp <- K[p, p, drop = FALSE]
    Lp <- t(if (inherits(K, "advector")) chol(Kp) else gp_chol_num(Kp))
    Lp[order(p), , drop = FALSE]
  }
  perm_dtheta <- function(z, theta, blk) {
    p <- ncgp_perm(blk)
    K <- gp_corr(theta, blk)
    Lp <- t(gp_chol_num(K[p, p, drop = FALSE]))
    out <- matrix(0, length(z), length(theta))
    out[, 1L] <- as.vector(Lp[order(p), , drop = FALSE] %*% z) * 0
    # the first column is u with sd L u = d f / d log sd, so z in this basis
    out[, 1L] <- z
    K0 <- K
    diag(K0) <- diag(K0) - gp_nugget
    for (k in seq_len(length(theta) - 1L)) {
      dK <- 0
      for (j in seq_along(blk[["aux_D2"]])) {
        if (isTRUE(blk[["gp_iso"]]) || j == k) {
          rho <- exp(theta[1L + (if (isTRUE(blk[["gp_iso"]])) 1L else j)])
          dK <- dK + blk[["aux_D2"]][[j]] / rho^2
        }
      }
      dK <- (K0 * dK)[p, p, drop = FALSE]
      M <- backsolve(Lp, t(backsolve(Lp, dK, upper.tri = FALSE)), upper.tri = FALSE)
      M[upper.tri(M)] <- 0
      diag(M) <- 0.5 * diag(M)
      out[, 1L + k] <- as.vector(M %*% z)
    }
    out
  }
  environment(perm_factor)$ncgp_perm <- ncgp_perm
  utils::assignInNamespace("gp_factor", perm_factor, ns = "frmtmb")
  utils::assignInNamespace("gp_field_dtheta", perm_dtheta, ns = "frmtmb")
})
