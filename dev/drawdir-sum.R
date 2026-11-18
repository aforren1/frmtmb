# Lane drawdir: summarise dev/drawdir-probe-*.tsv by design. Per design and
# per direction set and probe distance c: the largest and smallest
# per-direction excess (as the drop's ratio to the law's, e / (c^2 / 2)),
# over fits, split by the summed check's verdict.
# Usage: Rscript dev/drawdir-sum.R <tsv> [<tsv> ...]
args <- commandArgs(TRUE)
d <- do.call(rbind, lapply(args, function(f) {
  x <- utils::read.delim(f, stringsAsFactors = FALSE)
  x$file <- basename(f)
  x
}))
d$design <- sub(":[0-9]+$", "", d$id)
d$design <- sub("^(partR:[^:]+):.*$", "\\1", d$design)
d$design <- sub("^gglm:[0-9]+:[0-9]+$", "gglm", d$design)
d$design <- sub("^(gri|grs):.*$", "\\1", d$design)
d$design <- sub("^(wb|kd):", "", d$design)
d$design <- sub("^([A-Za-z_]+[0-9.]*):.*", "\\1", d$design)
cs <- as.numeric(strsplit(Sys.getenv("DRAWDIR_CS", "1.5,2,3"), ",")[[1]])
f3 <- function(x) formatC(x, digits = 3, format = "g")
for (kind in c("eig", "cor", "par")) {
  for (cc in cs) {
    tag <- paste0(kind, cc)
    r <- d[[paste0(tag, "_max")]] / (cc^2 / 2)
    s <- d[[paste0(tag, "_min")]] / (cc^2 / 2)
    cat(sprintf("\n== %s c %.1f: drop / law's drop - 1, max and min over directions\n",
                kind, cc))
    cat(sprintf("%-34s %4s %5s | %-26s | %-26s\n", "design", "fits", "fired",
                "unfired: max of max / min", "fired: min of max"))
    for (g in unique(d$design)) {
      j <- d$design == g
      fi <- j & d$fire
      un <- j & !d$fire
      cat(sprintf("%-34s %4d %5d | %8s %8s %8s | %8s\n", g, sum(j), sum(fi),
                  if (any(un)) f3(max(r[un], na.rm = TRUE)) else "-",
                  if (any(un)) f3(stats::quantile(r[un], 0.9, na.rm = TRUE))
                  else "-",
                  if (any(un)) f3(min(s[un], na.rm = TRUE)) else "-",
                  if (any(fi)) f3(min(r[fi], na.rm = TRUE)) else "-"))
    }
  }
}
cat("\n== two-sided reading of the 40 test draws (5644): q25 and q75\n")
for (g in unique(d$design)) {
  j <- d$design == g
  cat(sprintf("%-34s fits %4d | q75 max %7.3g  n(q75 > 1) %3d | q25 min %7.3g n(q25 < -0.5) %3d | fired %d\n",
              g, sum(j), max(d$q75[j]), sum(d$q75[j] > 1), min(d$q25[j]),
              sum(d$q25[j] < -0.5), sum(d$fire[j])))
}
