# Lane nlstop2, item 5079: where each Michaelis-Menten fit of the review's
# fresh rows sits against the poles of its body. For K the walls are at
# -xp; prints the fit's K, the band [-max xp, -min xp], the distance from K
# to the nearest wall and the median gap between walls.
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-polegeo.R
source("dev/nlstop2-common.R")
fresh_defs()
cat("ROWHEAD fam seed ll code wclass mclass Vm K band_lo band_hi dist gap\n")
for (fam in c("gaussian", "poisson", "bernoulli", "binomial", "gamma",
              "cumulative")) {
  for (seed in c(11, 22, 33, 44, 55)) {
    d <- gen(fam, "mm", seed * 1000 + match("mm", names(B)))
    r <- tt(frm(form(fam, "mm"), family = F[[fam]], data = d))
    if (iserr(r)) { cat("ROW", fam, seed, "ERR\n"); next }
    b <- unname(r$estimates$beta)
    nb <- names(r$estimates$beta)
    K <- b[grep("^K_", nb)]; Vm <- b[grep("^Vm_", nb)]
    w <- sort(unique(-d$xp))
    cat(sprintf("ROW %s %d %.6f %s %s %s %.5g %.6g %.4g %.4g %.3g %.3g\n",
                fam, seed, llv(r), r$opt$convergence, wsum(r), msum(r), Vm,
                K, min(w), max(w), min(abs(K - w)), median(diff(w))))
  }
}
cat("DONE\n")
