# Reviewer of the 0.70.0 consolidation: summarize dev/relrev070-interplay.R
# (release build r8 against lane moparity's build), per design and scale.
#   Rscript dev/relrev070-interplay-sum.R
L <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev070-log/inter"
a <- utils::read.delim(file.path(L, "r8.tsv"), stringsAsFactors = FALSE)
b <- utils::read.delim(file.path(L, "lane.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(a$seed, b$seed), identical(a$design, b$design),
          identical(a$scale, b$scale))
d <- a$ll - b$ll
a$warn[is.na(a$warn)] <- ""
b$warn[is.na(b$warn)] <- ""
key <- paste(a$design, a$scale)
for (k in unique(key)) {
  i <- key == k
  cat(sprintf(paste0("%-14s n %2d | r8 - lane logLik: identical %2d, ",
                     "> +1e-6 %2d, < -1e-6 %2d (min %.3g, max %.3g) | ",
                     "code!=0 r8 %d lane %d | warned r8 %d lane %d | ",
                     "NaN-SE fits r8 %d lane %d | polished r8 %d | ",
                     "mo runs r8 %d lane %d | kappa max %.3g | spread %s\n"),
              k, sum(i), sum(d[i] == 0), sum(d[i] > 1e-6), sum(d[i] < -1e-6),
              min(d[i]), max(d[i]), sum(a$code[i] != 0), sum(b$code[i] != 0),
              sum(a$nwarn[i] > 0), sum(b$nwarn[i] > 0),
              sum(a$nse_nan[i] > 0), sum(b$nse_nan[i] > 0),
              sum(a$polish[i] > 0), sum(a$mo_runs[i]), sum(b$mo_runs[i]),
              max(a$kappa[i]), paste(unique(a$spread_mo[i]), collapse = ";")))
}
cat("\nRows where the warning text differs between builds:\n")
w <- which(a$warn != b$warn)
for (j in w) {
  cat(sprintf("  seed %d %s x%g | dll %.3g | r8: %s\n      lane: %s\n",
              a$seed[j], a$design[j], a$scale[j], d[j],
              substr(a$warn[j], 1, 160), substr(b$warn[j], 1, 160)))
}
cat("\nRows where r8 is worse than lane by more than 1e-6:\n")
for (j in which(d < -1e-6)) {
  cat(sprintf("  seed %d %s x%g | r8 %.8f code %d | lane %.8f code %d\n",
              a$seed[j], a$design[j], a$scale[j], a$ll[j], a$code[j],
              b$ll[j], b$code[j]))
}
cat("\nScale invariance within a build (ll at x1000 minus ll at x1, adjusted):\n")
for (bb in list(list("r8", a), list("lane", b))) {
  X <- bb[[2]]
  for (ds in unique(X$design[X$scale == 1000])) {
    u <- X[X$design == ds & X$scale == 1, ]
    v <- X[X$design == ds & X$scale == 1000, ]
    dd <- v$ll - u$ll
    cat(sprintf("  %-5s %-8s max |d| %.3g, > 1e-6 on %d of %d\n", bb[[1]], ds,
                max(abs(dd)), sum(abs(dd) > 1e-6), length(dd)))
  }
}
