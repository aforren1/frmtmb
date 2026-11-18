# Review of lane nbrobust: read dev/reviews/nbrobust-review-formula.tsv
t <- read.delim("dev/reviews/nbrobust-review-formula.tsv")
for (a in c("scalar", "robust")) {
  s <- t[t$arm == a, ]
  cat("worst relative value error,", a, "\n")
  print(head(s[order(-s$relv), c("y", "lmu", "lsh", "v", "vref", "relv",
                                 "absv")], 6))
}
s <- t[t$arm == "scalar", ]
s$es <- abs(s$vstats - s$vref)
cat("stats::dnbinom, worst absolute error\n")
print(head(s[order(-s$es), c("y", "lmu", "lsh", "vstats", "vref", "v",
                             "es")], 12))
cat("rows with stats err > 1e-8:", sum(s$es > 1e-8), "of", nrow(s), "\n")
print(table(lsh = s$lsh[s$es > 1e-8]))
print(table(y = s$y[s$es > 1e-8]))
cat("nb_lpmf rel err > 1e-11 rows:\n")
print(s[s$relv > 1e-11, c("y", "lmu", "lsh", "v", "vref", "relv")])
