# Review of lane nbrobust: the max relative difference per regime and
# family in dev/reviews/nbrobust-review-parity-<arm>.log
for (a in c("base", "lib")) {
  t <- readLines(paste0("dev/reviews/nbrobust-review-parity-", a, ".log"))
  t <- t[grepl("max rel", t)]
  r <- as.numeric(sub(".*max rel +([^ ]+) .*", "\\1", t))
  fam <- sub("^.*y +[0-9]+ ([a-z_0-9]+) .*$", "\\1", t)
  reg <- sub(" .*", "", t)
  cat("==", a, "\n")
  print(signif(tapply(r, list(fam, reg), max), 3))
}
