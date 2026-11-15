# Summary of dev/predraw2-proto.R's TSVs by design.
# Usage: Rscript dev/predraw2-proto-sum.R <tsv>...
fs <- commandArgs(TRUE)
d <- do.call(rbind, lapply(fs, function(f) {
  x <- utils::read.delim(f, stringsAsFactors = FALSE)
  x$file <- basename(f)
  x
}))
d$kind <- sub(":.*", "", d$id)
d$kind[d$kind == "sweep"] <- paste0("sweep_", sub(".*:", "", d$id[d$kind == "sweep"]))
qs <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return("-")
  sprintf("%.3f/%.3f/%.3f", min(v), stats::median(v), max(v))
}
cat(sprintf("%-16s %4s %-24s %-24s %-24s %-24s %-24s %s\n", "design", "n",
            "med20 min/med/max", "max20", "med1000", "k", "ess", "n2_20>0"))
for (k in unique(d$kind)) {
  x <- d[d$kind == k, ]
  cat(sprintf("%-16s %4d %-24s %-24s %-24s %-24s %-24s %d\n", k, nrow(x),
              qs(x$med20), qs(x$max20), qs(x$med), qs(x$k), qs(x$ess),
              sum(x$n2_20 > 0, na.rm = TRUE)))
}
