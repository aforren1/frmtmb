# Merge check 9: compare the confint_varcorr() rows of two dumps of
# mc9-confint-dump.R, per case group.
#   Rscript dev/rel083-kt/arspat2/mc9-cmp.R <tagA> <tagB>
a <- commandArgs(TRUE)
rd <- function(t) readRDS(sprintf("dev/rel083-log/kt/arspat2/mc9-%s.rds", t))
A <- rd(a[1]); B <- rd(a[2])
grp <- function(x) sub(" (k|seed) .*$", "", sub("^ou times.*", "ou", x))
common <- intersect(unique(A$case), unique(B$case))
cat(a[1], "cases", length(unique(A$case)), "|", a[2], "cases",
    length(unique(B$case)), "| common", length(common), "\n")
for (g in unique(grp(common))) {
  cs <- common[grp(common) == g]
  nid <- 0L; mx <- 0; nadiff <- 0L; ncase <- length(cs)
  for (c in cs) {
    x <- A[A$case == c, ]; y <- B[B$case == c, ]
    if (nrow(x) != nrow(y) || !identical(x$term, y$term)) {
      cat("  case", c, ": rows or terms differ\n"); next
    }
    m <- as.matrix(x[, c("estimate", "lwr", "upr")])
    n <- as.matrix(y[, c("estimate", "lwr", "upr")])
    if (identical(m, n)) nid <- nid + 1L
    nadiff <- nadiff + sum(is.na(m) != is.na(n))
    d <- abs(m - n); d <- d[is.finite(d)]
    if (length(d)) mx <- max(mx, d)
  }
  cat(sprintf("%-24s cases %3d identical %3d  NA pattern differs in %d cells  max |diff| %.3g\n",
              g, ncase, nid, nadiff, mx))
}
# ci082's expectations on the fixed correlations, per dump
for (t in a) {
  X <- rd(t)
  fx <- X[grepl("^xc[sh] k", X$case) & X$type == "cor", ]
  ok <- identical(fx$lwr, fx$estimate) && identical(fx$upr, fx$estimate)
  cat(sprintf("%s: xcs/xch fixed cor rows %d, lwr == upr == estimate on %d, NA bounds on %d -> %s\n",
              t, nrow(fx), sum(fx$lwr == fx$estimate & fx$upr == fx$estimate,
                               na.rm = TRUE), sum(is.na(fx$lwr)),
              if (ok) "ci082 expectation holds" else "ci082 expectation FAILS"))
  sr <- X[grepl("^sar", X$case), ]
  if (nrow(sr)) {
    cat(sprintf("%s: sar rows %d, types %s, coverage of 0.9 %.3f, upr at 1 or above %d, NA %d\n",
                t, nrow(sr), paste(unique(sr$type), collapse = ","),
                mean(sr$lwr <= 0.9 & 0.9 <= sr$upr, na.rm = TRUE),
                sum(sr$upr >= 1, na.rm = TRUE), sum(is.na(sr$lwr))))
  }
}
