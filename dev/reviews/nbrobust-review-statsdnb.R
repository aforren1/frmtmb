# Review of lane nbrobust: stats::dnbinom(mu =) against the mpfr
# reference at moderate means, from dev/reviews/nbrobust-review-formula.tsv
t <- read.delim("dev/reviews/nbrobust-review-formula.tsv")
s <- t[t$arm == "scalar" & t$lmu <= 1.5 & t$lmu >= -1, ]
s$es <- abs(s$vstats - s$vref) / abs(s$vref)
s$en <- abs(s$v - s$vref) / abs(s$vref)
cat("lmu in [-1, 1.5]: max rel error by log shape, stats::dnbinom and nb_lpmf\n")
print(signif(cbind(stats = tapply(s$es, s$lsh, max),
                   nb_lpmf = tapply(s$en, s$lsh, max)), 3))
