# rel084-kt copy of dev/reviews/alapse-review-qglob.R: sources dev/rel084-kt/alapse-review-common.R (rellib-r23). No other change.
# Review of lane alapse, check 1 (the filed MEDIUM): an estimated
# quantile proved at a local maximum below the global one. For design B
# of dev/alapse-rq.R (100 rows, y = 1 + 0.5 x + N(0, 1), y ~ x, quantile
# estimated), fresh seeds, the exact global maximum over the quantile
# process (every rq() solution; the profile in p is smooth on each
# interval) against the fit's logLik; and on design E (quantile ~ z)
# multi-start alternating maximization.
#   Rscript dev/reviews/alapse-review-qglob.R <base|lane> <seeds>
a <- commandArgs(TRUE)
lib_arg <- a[1]
seeds <- eval(parse(text = a[2]))
source("dev/rel084-kt/alapse-review-common.R")

exact_const <- function(y, X) {
  n <- length(y)
  pr <- suppressWarnings(quantreg::rq(y ~ X - 1, tau = -1))
  sol <- pr$sol
  taus <- sol[1, ]
  B <- sol[-(1:3), , drop = FALSE]
  best <- list(ll = -Inf)
  for (j in seq_len(ncol(B))) {
    lo <- max(taus[j], 1e-9)
    hi <- min(if (j < ncol(B)) taus[j + 1] else 1, 1 - 1e-9)
    if (hi <= lo) next
    r <- y - drop(X %*% B[, j])
    A <- sum(pmax(r, 0)); Bn <- sum(pmax(-r, 0))
    prof <- function(p) n * log(p) + n * log1p(-p) -
      n * log((p * A + (1 - p) * Bn) / n) - n
    o <- optimize(prof, c(lo, hi), maximum = TRUE, tol = 1e-12)
    for (pp in c(o$maximum, lo, hi)) {
      v <- prof(pp)
      if (is.finite(v) && v > best$ll) best <- list(ll = v, p = pp)
    }
  }
  best
}

cat("design B\n")
res <- NULL
for (s in seeds) {
  set.seed(2000 + s)
  d <- data.frame(x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  r <- cap(frm(y ~ x, family = asym_laplace(), data = d))
  f <- r$v
  p <- ns$eval_dpars(f)[[1]]$quantile[1]
  ex <- exact_const(d$y, cbind(1, d$x))
  ll <- as.numeric(logLik(f))
  conv <- any(r$wc %in% c("frmtmb_conv_code", "frmtmb_conv_gradient"))
  cat(sprintf("B %3d code %d cert %-5s conv-warn %-5s ll %.6f global %.6f gap %.3e p %.4f p* %.4f\n",
              s, f$opt$convergence, isTRUE(f$opt$ald_cert$ok), conv, ll,
              ex$ll, ex$ll - ll, p, ex$p))
  res <- rbind(res, data.frame(s = s, gap = ex$ll - ll, dp = abs(p - ex$p),
                               cert = isTRUE(f$opt$ald_cert$ok),
                               conv = conv))
}
cat(sprintf("B summary: n %d | gap > 1e-6: %d | proved and gap > 1e-6: %d | proved, gap > 1e-6, no convergence warning: %d | max gap %.3g | max dp %.3g\n",
            nrow(res), sum(res$gap > 1e-6), sum(res$cert & res$gap > 1e-6),
            sum(res$cert & res$gap > 1e-6 & !res$conv), max(res$gap),
            max(res$dp)))
cat("DONE\n")
