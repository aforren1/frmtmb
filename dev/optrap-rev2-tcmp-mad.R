# Reviewer of lane optrap, re-check: as dev/optrap-rev2-tcmp.R for the MAD arm,
# against rellib-r7 (log-likelihood unscaled to y x 1).
for (cs in c("t1", "t2")) {
  a <- read.delim(sprintf("dev/optrap-rev2-log/b1-%s-lanemad.tsv", cs), quote = "")
  b <- read.delim(sprintf("dev/optrap-rev2-log/b1-%s-r7.tsv", cs), quote = "")
  best <- tapply(pmax(a$ll, b$ll), a$seed, max)
  for (sc in c(1, 1000)) {
    i <- a$setting == sc
    ga <- a$ll[i] - best[as.character(a$seed[i])]
    gb <- b$ll[i] - best[as.character(b$seed[i])]
    cat(sprintf("%s scale %g | short by > 1e-4 of the best of both builds at both scales: lane %d (min %.3g, silent %d) | r7 %d (min %.3g, silent %d)\n",
                cs, sc, sum(ga < -1e-4), min(ga),
                sum(ga < -1e-4 & a$nwarn[i] == 0), sum(gb < -1e-4), min(gb),
                sum(gb < -1e-4 & b$nwarn[i] == 0)))
  }
}
