# Reviewer of lane illse: fitted() SEs on raw polynomials with (1 | g)
# whose sd is at zero (route skipped), against lm() (the sd at zero is
# the lm fit; ML scaling).
#   REV_ARM=lane|base Rscript dev/illse-rev-edge2.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
for (deg in c(6, 7, 8)) for (seed in 1:4) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
  fo0 <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
  r <- cap(frm(fo, data = d))
  f <- r$v
  sd <- exp(f$opt$par[names(f$opt$par) == "theta"])
  if (sd > 1e-4) { cat(sprintf("deg %d s%d sd %.2e not at zero\n", deg, seed, sd)); next }
  lmf <- stats::lm(fo0, data = d)
  nd <- data.frame(x = c(1.1, 1.5, 1.9), g = d$g[1])
  for (k in seq_len(deg)) nd[[xn[k]]] <- nd$x^k
  ref <- stats::predict(lmf, nd, se.fit = TRUE)$se.fit * sqrt((n - deg - 1) / n)
  fr <- cap(fitted(f, newdata = nd, re_formula = NA))
  se <- if (inherits(fr$v, "rev_err")) NA else fr$v[, "Est.Error"]
  l <- ns$sdr_of(f)$se_lost
  cat(sprintf("deg %d s%d | sd %.1e | route %s | lost %s | fit warn %s | fitted SE / lm %s | fitted warn %s | ll short %.2e\n",
              deg, seed, sd, !is.null(ns$sdr_of(f)$cov_fixed_root),
              if (length(l)) paste(paste0(names(l), ":", l), collapse = " ") else "-",
              wtag(r$w), ratio_str(se / ref), wtag(fr$w),
              as.numeric(logLik(lmf)) - as.numeric(-f$opt$objective - 0)))
}
