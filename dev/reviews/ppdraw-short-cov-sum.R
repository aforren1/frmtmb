# Summary of ppdraw-short-cov.R: lane against 0.81.0 against the plug-in,
# on the fits the side law changes and on all fits.
# Usage: Rscript ppdraw-short-cov-sum.R <lane.tsv> <base.tsv>
args <- commandArgs(TRUE)
ln <- utils::read.delim(args[1])
bs <- utils::read.delim(args[2])
key <- c("case", "par", "seed")
L <- ln[ln$pe, ]
B <- bs[bs$pe, ]
P <- bs[!bs$pe, ]
m <- merge(merge(L, B, by = key, suffixes = c(".lane", ".base")),
           P[, c(key, "cov90", "w90", "cov98", "w98", "p0", "est")],
           by = key)
for (cs in c("R", "S")) for (pr in unique(m$par[m$case == cs])) {
  s <- m[m$case == cs & m$par == pr, ]
  for (sub in c("changed", "all")) {
    t <- if (sub == "changed") s[s$changed.lane, ] else s
    if (!nrow(t)) next
    cat(sprintf(paste("%s par %.2f %-7s n %3d | cov90 lane %.3f base %.3f",
                      "plug %.3f | cov98 lane %.3f base %.3f plug %.3f |",
                      "median w98 lane %.3f base %.3f plug %.3f\n"),
                cs, pr, sub, nrow(t), mean(t$cov90.lane), mean(t$cov90.base),
                mean(t$cov90), mean(t$cov98.lane), mean(t$cov98.base),
                mean(t$cov98), stats::median(t$w98.lane),
                stats::median(t$w98.base), stats::median(t$w98)))
  }
  t <- s[s$changed.lane, ]
  if (nrow(t)) {
    cat(sprintf("   changed fits where lane cov98 < base cov98: %d of %d;",
                sum(t$cov98.lane < t$cov98.base), nrow(t)),
        sprintf("lane cov98 < 0.98: %d; base cov98 < 0.98: %d\n",
                sum(t$cov98.lane < 0.98), sum(t$cov98.base < 0.98)))
  }
}
z <- m[m$case == "Z", ]
for (sub in c("changed", "all")) {
  t <- if (sub == "changed") z[z$changed.lane, ] else z
  cat(sprintf(paste("Z %-7s n %2d | mean |P0 - true| lane %.4f base %.4f",
                    "plug %.4f | mean |mean - true| lane %.4f base %.4f",
                    "plug %.4f\n"),
              sub, nrow(t), mean(abs(t$p0.lane - t$p0true.lane)),
              mean(abs(t$p0.base - t$p0true.lane)),
              mean(abs(t$p0 - t$p0true.lane)),
              mean(abs(t$est.lane - t$truth.lane)),
              mean(abs(t$est.base - t$truth.lane)),
              mean(abs(t$est - t$truth.lane))))
}
