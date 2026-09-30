# SUPERSEDED in punch 2: written against the punch-1
# log_ibeta_cf(x, a, b, N, lx, l1mx), whose signature punch 2 changed,
# so it no longer runs on the lane build. Its output is kept.
# Punch 1, B2: how far past m does each fraction stay exact? The direct
# fraction for I_x(a, b) and the complement for 1 - I_{1-x}(b, a), the
# latter reading log(x) and log(1 - x) exactly, at x = m + k sd, against
# stats::pbeta(log.p = TRUE) (relative, floored at one). 50 steps.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cf <- frmtmb:::log_ibeta_cf
ks <- c(-2, -1, -0.5, 0, 0.5, 1, 1.25, 1.5, 2)
rows <- list()
for (a in c(0.01, 0.3, 4, 60, 400)) for (b in c(0.5, 7, 300, 7e3, 8e5)) {
  s <- a + b; m <- (a + 1) / (s + 2); sd <- sqrt(a * b / (s * s * (s + 1)))
  for (k in ks) {
    x <- m + k * sd
    if (x <= 0 || x >= 0.5) next
    r <- pbeta(x, a, b, log.p = TRUE)
    if (!is.finite(r) || r < -600) next
    d <- cf(x, a, b, 50L)
    cmp <- log1p(-exp(cf(1 - x, b, a, 50L, lx = log1p(-x), l1mx = log(x))))
    rows[[length(rows) + 1]] <- data.frame(a, b, k, m = signif(m, 3),
      direct = abs(d - r) / max(1, abs(r)), comp = abs(cmp - r) / max(1, abs(r)))
  }
}
res <- do.call(rbind, rows)
options(width = 150)
res$direct <- signif(res$direct, 2); res$comp <- signif(res$comp, 2)
print(reshape(res[, c("a", "b", "k", "direct")], idvar = c("a", "b"), timevar = "k", direction = "wide"))
print(reshape(res[, c("a", "b", "k", "comp")], idvar = c("a", "b"), timevar = "k", direction = "wide"))
