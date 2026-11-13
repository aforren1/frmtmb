# Short check of lane ci082's fix round: how far a FIXED correlation read
# back through two matrix inverses moves between varcorr_fixed_rows()'s
# probe points, by size of the fixed AR(1) correlation matrix, with no
# fit (theta = log sds spread from -3 to 7). Lane library only.
#   Rscript dev/reviews/ci082-short-invprobe.R
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
# AR(1) R (rr < 1), then a gaussian kernel exp(-(i - j)^2 / (2 l^2)) with
# l = -rr (rr < 0): no correlation past the 0.9999 clamp, far worse
# conditioned
cases <- rbind(expand.grid(dd = c(6L, 10L, 20L, 40L),
                           rr = c(0.995, 0.999, 0.9995)),
               expand.grid(dd = c(6L, 10L, 20L), rr = c(-1, -2, -3)))
for (ci in seq_len(nrow(cases))) {
  dd <- cases$dd[ci]
  rr <- cases$rr[ci]
  lag <- abs(outer(seq_len(dd), seq_len(dd), "-"))
  R <- if (rr > 0) rr^lag else exp(-lag^2 / (2 * rr^2))
  nm <- sprintf("xi%d_%d", dd, ci)
  frmtmb_register_covstruct(
    nm,
    lpdf = function(b, theta, blk) 0,
    vcov = local({
      R <- R
      function(theta, blk) {
        n <- blk[["dim"]]
        s <- exp(theta[seq_len(n)])
        Q <- diag(1 / s, n) %*% solve(R) %*% diag(1 / s, n)
        V <- solve(Q)
        dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
        V
      }
    }),
    npar = dd, sd_idx = seq_len(dd))
  bk <- list(covstruct = nm, dim = dd, cnms = sprintf("t%d", seq_len(dd)),
             theta_idx = seq_len(dd))
  t0 <- seq(-3, 7, length.out = dd)
  fx <- ns$varcorr_fixed_rows(bk, t0)
  nsd <- dd
  # the largest relative move of a correlation row over the +/-1 points,
  # and the largest error of the read-back correlation against R
  e <- ns$cs_entry(nm)
  rows <- function(tt) {
    C <- stats::cov2cor(e$vcov(tt, bk))
    atanh(C[lower.tri(C)])
  }
  g0 <- tryCatch(rows(t0), error = function(err) NULL)
  if (is.null(g0)) {
    cat(sprintf("RESULT d %d rr %g: vcov errors at t0; fixed is %s\n", dd,
                rr, if (is.null(fx)) "NULL (nothing fixed)" else "set"))
    next
  }
  mv <- 0
  for (i in seq_along(t0)) for (s in c(-1, 1)) {
    tt <- t0
    tt[i] <- tt[i] + s
    gp <- rows(tt)
    mv <- max(mv, abs(gp - g0) / pmax(abs(g0), abs(gp)), na.rm = TRUE)
  }
  err <- max(abs(tanh(g0) - R[lower.tri(R)]))
  cat(sprintf(paste("RESULT d %d %s kappa %.3g max cor %.6f cor rows %d",
                    "fixed %d max-rel-move %.3g abs-err-vs-R %.3g\n"), dd,
              if (rr > 0) sprintf("ar1 %.4f", rr) else
                sprintf("gauss l %g", -rr),
              kappa(R, exact = TRUE), max(R[lower.tri(R)]),
              length(fx) - nsd, sum(fx[-seq_len(nsd)]), mv, err))
}
