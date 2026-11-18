# Lane drawdir, backlog 5644: a two-sided reading of the summed check's 40
# test draws (a quantile on each side, not the median), from the
# quantiles dev/drawdir-probe.R records. Per design: how many fits it
# fires on at q75 > 1 and q25 < -0.5 (the item's S0.1 seed 14 reads 4.03
# and -0.52 over 1000 draws), and at q75 > 1 alone.
# Usage: Rscript dev/drawdir-twoside.R <probe tsv> [...]
d <- do.call(rbind, lapply(commandArgs(TRUE), function(f) {
  x <- utils::read.delim(f, stringsAsFactors = FALSE)
  x[c("id", "fire", "q25", "q75")]
}))
d$design <- sub(":[0-9]+$", "", d$id)
d$design <- sub("^(gglm:[0-9]+):[0-9]+$", "\\1", d$design)
d$design <- sub("^(gri|grs):.*$", "\\1", d$design)
d$design <- sub("^([A-Za-z_]+[0-9.]*):.*", "\\1", d$design)
d$both <- d$q75 > 1 & d$q25 < -0.5
cat(sprintf("%-22s %5s %8s %10s %12s\n", "design", "fits", "summed",
            "q75>1", "q75>1&q25<-.5"))
for (g in unique(d$design)) {
  j <- d$design == g
  cat(sprintf("%-22s %5d %8d %10d %12d\n", g, sum(j), sum(d$fire[j]),
              sum(d$q75[j] > 1), sum(d$both[j])))
}
cat(sprintf("all %d fits: q75>1 %d, both %d; S0.1:14 both %s\n", nrow(d),
            sum(d$q75 > 1), sum(d$both),
            paste(d$both[d$id == "S0.1:14"], collapse = ",")))
