# Short check of lane squant fix round 1, item 3, independent reference:
# squant-short-quant.R's round trip read the package's own elf_cdf_log()
# (which the Newton steps solve on) and pbeta(), which is -Inf past
# |t| = 745. Here pr' at the returned q is from numerical integration of
# an independently written density, in log space (squant-short-tape.R's
# method, with breakpoints out to 2^25 for the slow tails).
# Run (worktree root):
#   Rscript dev/reviews/squant-short-quant2.R > dev/reviews/squant-short-quant2.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
eq <- asNamespace("frmtmb")$elf_quantile
lpdf_ref <- function(y, mu, sigma, p, lam) {
  e <- (y - mu) / sigma
  z <- e / lam
  sp <- ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
  (1 - p) * e - lam * sp - log(lam * sigma) - lbeta(lam * (1 - p), lam * p)
}
tail_int <- function(q, s, mu, sigma, p, lam) {
  l0 <- lpdf_ref(q, mu, sigma, p, lam)
  br <- c(0, sigma * 2^seq(-24, 25, by = 1))
  tot <- 0
  for (k in seq_len(length(br) - 1L)) {
    v <- stats::integrate(function(v) {
      exp(lpdf_ref(q + s * v, mu, sigma, p, lam) - l0)
    }, br[k], br[k + 1L], rel.tol = 1e-12, subdivisions = 1000L)$value
    tot <- tot + v
    if (tot > 0 && v < 1e-17 * tot && br[k + 1L] > 50 * sigma / min(p, 1 - p))
      break
  }
  l0 + log(tot)
}
ref_side <- function(q, mu, sigma, p, lam, up) {
  mode <- mu + lam * sigma * stats::qlogis(1 - p)
  if (!up) {
    if (q <= mode) tail_int(q, -1, mu, sigma, p, lam) else
      log1p(-exp(tail_int(q, 1, mu, sigma, p, lam)))
  } else {
    if (q > mode) tail_int(q, 1, mu, sigma, p, lam) else
      log1p(-exp(tail_int(q, -1, mu, sigma, p, lam)))
  }
}
prs <- c(1e-300, 1e-100, 1e-20, 1e-8, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-3,
         1 - 1e-8, 1 - 1e-15)
worst <- 0
n <- 0
nbad <- 0
for (p in c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)) {
  for (lam in c(0.001, 0.01, 0.05, 0.25, 1, 4)) {
    for (ms in list(c(0, 1), c(3, 2))) {
      q <- eq(prs, ms[1], ms[2], p, lam)
      up <- prs > 0.5
      goal <- ifelse(up, log1p(-prs), log(prs))
      got <- mapply(ref_side, q, up, MoreArgs = list(mu = ms[1],
                                                       sigma = ms[2],
                                                       p = p, lam = lam))
      # relative error of the tail probability, and of the log where
      # the probability is near 1 (pr 0.5)
      err <- abs(got - goal)
      worst <- max(worst, err)
      n <- n + length(err)
      nbad <- nbad + sum(!is.finite(err) | err > 1e-8)
      cat(sprintf("p %.3f lam %.3f mu %g sigma %g: max |log tail - goal| %.2g at pr %s%s\n",
                  p, lam, ms[1], ms[2], max(err),
                  format(prs[which.max(err)], digits = 3),
                  if (any(!is.finite(err) | err > 1e-8))
                    paste0(" | BAD at pr ", paste(format(prs[!is.finite(err) | err > 1e-8], digits = 3), collapse = ","))
                  else ""))
    }
  }
}
cat(sprintf("SUMMARY: worst |log tail(q) - log pr| %.3g over %d points; above 1e-8 or nonfinite: %d\n",
            worst, n, nbad))
cat("DONE\n")
