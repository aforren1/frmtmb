# Reviewer: a cumulative("probit") row with no cs() whose two thresholds
# sit a few ulps apart near 0.67449, where pnorm() is not monotone in
# floating point (dev/relrev-mono.R): ord_sim() on rellib-r5 and r6.
for (lib in c("C:/Users/adf44/source/r/rellib-r5", "C:/Users/adf44/source/r/rellib-r6")) {
  loadNamespace("frmtmb", lib.loc = c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  fam <- frmtmb::cumulative("probit")
  a <- 0.67448975
  x <- a + (seq_len(200000) - 100000) * a * .Machine$double.eps
  hit <- NULL
  for (i in which(diff(pnorm(x)) < 0)) {
    for (g in 1:4) {
      raw <- c(x[i], log(g * 2.220446e-16), log(1))
      tau <- frmtmb:::ord_threshold_values(fam, raw)
      if (tau[2] > tau[1] && pnorm(tau[2]) < pnorm(tau[1])) { hit <- raw; break }
    }
    if (!is.null(hit)) break
  }
  tau <- frmtmb:::ord_threshold_values(fam, hit)
  dp <- list(mu = 0, disc = 1)
  set.seed(1)
  s <- vapply(1:200, function(j) as.numeric(fam$sim(dp, list(), 1L, list(tau_raw = hit))), 0)
  lp <- fam$lpdf(c(1, 2, 3, 4), list(mu = rep(0, 4), disc = rep(1, 4)), list(), list(tau_raw = hit))
  cat(format(packageVersion("frmtmb", lib.loc = lib)), ": tau2 - tau1", tau[2] - tau[1],
      " P2", pnorm(tau[2]) - pnorm(tau[1]), " NA draws", sum(is.na(s)), "of 200;",
      " lpdf y=1..4:", format(lp, digits = 4), "\n")
  unloadNamespace("frmtmb")
}
