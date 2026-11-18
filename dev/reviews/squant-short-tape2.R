# Short check of lane squant fix round 1, item 1 follow-up to
# squant-short-tape.txt: (a) relative error of the taped d lF / d mu and
# d lS / d mu where the analytic value is between 1e-30 and 1e-12 of the
# slope scale max(p, 1 - p) (the band the first script left as 'tiny');
# (b) the 11 taped d/dt values across t = -690 and 690 where that script
# reported a relative jump of 1, with the analytic value beside each.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-tape2.R > dev/reviews/squant-short-tape2.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
lpdf_ref <- function(y, mu, sigma, p, lam) {
  e <- (y - mu) / sigma
  z <- e / lam
  sp <- ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
  (1 - p) * e - lam * sp - log(lam * sigma) - lbeta(lam * (1 - p), lam * p)
}
# the small tail on each side has the closed form of the leading term
# times a series; use integration as in squant-short-tape.R
tail_int <- function(q, s, p, lam) {
  l0 <- lpdf_ref(q, 0, 1, p, lam)
  br <- c(0, lam * 2^seq(-12, 22, by = 1))
  tot <- 0
  for (k in seq_len(length(br) - 1L)) {
    v <- stats::integrate(function(v) exp(lpdf_ref(q + s * v, 0, 1, p, lam) - l0),
                          br[k], br[k + 1L], rel.tol = 1e-13,
                          subdivisions = 500L)$value
    tot <- tot + v
    if (v < 1e-17 * tot) break
  }
  l0 + log(tot)
}
# on the large side, F (or S) is 1 - (small tail) and d log F = -f / F
# with F = 1 to double precision, so the analytic value is -f (or f)
an_d <- function(q, p, lam, side) {
  lf <- lpdf_ref(q, 0, 1, p, lam)
  mode <- lam * stats::qlogis(1 - p)
  if (side == "lF") {
    lF <- if (q <= mode) tail_int(q, -1, p, lam) else
      log1p(-exp(tail_int(q, 1, p, lam)))
    -exp(lf - lF)
  } else {
    lS <- if (q > mode) tail_int(q, 1, p, lam) else
      log1p(-exp(tail_int(q, -1, p, lam)))
    exp(lf - lS)
  }
}
q_of_t <- function(t, lam) {
  q <- lam * t
  if (q / lam == t) return(q)
  for (k in 1:8) for (s in c(-1, 1)) {
    qq <- q + s * k * .Machine$double.eps * abs(q)
    if (qq / lam == t) return(qq)
  }
  q
}
cat("== (a) band 1e-30 to 1e-12 of the slope scale ==\n")
tg <- sort(unique(c(seq(-689, 689, by = 1), -689.999999, 689.999999)))
for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) for (lam in c(0.05, 0.25, 1, 4)) {
  fam <- quantile_elf(smoothing = lam)
  scale <- max(p, 1 - p)
  for (side in c("lF", "lS")) {
    # the band lies on the large side only: pick t there
    tt <- if (side == "lF") tg[tg > 0] else tg[tg < 0]
    q <- vapply(tt, q_of_t, 0, lam = lam)
    an <- vapply(q, an_d, 0, p = p, lam = lam, side = side)
    band <- abs(an) > 1e-30 * scale & abs(an) <= 1e-12 * scale
    if (!any(band)) {
      cat(sprintf("p %.2f lam %.2f %s: no point in the band\n", p, lam, side))
      next
    }
    q <- q[band]
    an <- an[band]
    tt <- tt[band]
    tp <- RTMB::MakeTape(function(m) {
      fam$cdf_log(q, list(mu = m + 0 * q, sigma = 1, quantile = p),
                  list())[[side]]
    }, 0)
    g <- as.numeric(tp$jacobian(0))
    re <- abs(g / an - 1)
    cat(sprintf(paste0("p %.2f lam %.2f %s: %d points, t %g to %g, worst rel %.2g at t %g",
                       " (an %.3g, tape %.3g), n rel > 1e-6: %d\n"),
                p, lam, side, length(q), min(tt), max(tt), max(re),
                tt[which.max(re)], an[which.max(re)], g[which.max(re)],
                sum(re > 1e-6)))
  }
}
cat("\n== (b) the 11 points across +-690 ==\n")
for (cs in list(list(0.5, 0.25, "lF", 690), list(0.5, 0.25, "lS", -690),
                list(0.1, 1, "lF", 690), list(0.01, 0.25, "lS", -690),
                list(0.99, 1, "lF", 690))) {
  p <- cs[[1]]
  lam <- cs[[2]]
  side <- cs[[3]]
  t0 <- cs[[4]]
  fam <- quantile_elf(smoothing = lam)
  tt <- t0 + (-5:5) * 1e-6
  q <- vapply(tt, q_of_t, 0, lam = lam)
  tp <- RTMB::MakeTape(function(m) {
    fam$cdf_log(q, list(mu = m + 0 * q, sigma = 1, quantile = p),
                list())[[side]]
  }, 0)
  g <- as.numeric(tp$jacobian(0))
  an <- vapply(q, an_d, 0, p = p, lam = lam, side = side)
  cat(sprintf("p %.2f lam %.2f %s t0 %g: d/dmu tape %s\n  analytic %s\n",
              p, lam, side, t0, paste(sprintf("%.3g", g), collapse = " "),
              paste(sprintf("%.3g", an), collapse = " ")))
}
cat("DONE\n")
