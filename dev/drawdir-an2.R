# Lane drawdir: candidate per-direction statistics from
# dev/drawdir-probe2-*.tsv, per fit and probe distance c:
#   S1  max over parameters and sides of eA (the direction's excess)
#   S2  max of eA - max(eB, 0) on the same side (the direction's excess
#       beyond its own coordinate's slice)
# Printed by design: max over unfired fits; and the per-fit values of the
# fits named by a regex.
# Usage: Rscript dev/drawdir-an2.R <regex of fits to list> <tsv> [...]
a <- commandArgs(TRUE)
d <- do.call(rbind, lapply(a[-1], utils::read.delim, stringsAsFactors = FALSE))
d$S1 <- pmax(d$eAp, d$eAm)
d$S2 <- pmax(d$eAp - pmax(d$eBp, 0), d$eAm - pmax(d$eBm, 0))
f <- do.call(rbind, lapply(split(d, list(d$id, d$c), drop = TRUE), function(x) {
  data.frame(id = x$id[1], c = x$c[1], fire = x$fire[1], excess = x$excess[1],
             S1 = max(x$S1), S1par = x$par[which.max(x$S1)],
             S2 = max(x$S2), S2par = x$par[which.max(x$S2)])
}))
f$design <- sub(":[0-9]+$", "", f$id)
f$design <- sub("^(gri|grs):.*$", "\\1", f$design)
f$design <- sub("^(partR:[^:]+:[0-9]+):.*$", "\\1", f$design)
for (cc in sort(unique(f$c))) {
  cat(sprintf("\n== c %.1f: max over unfired fits | min over fired fits\n", cc))
  for (g in unique(f$design)) {
    x <- f[f$design == g & f$c == cc, ]
    u <- x[!x$fire, ]
    v <- x[x$fire, ]
    cat(sprintf("%-36s n %3d | S1 %8.3g S2 %8.3g | fired %2d S1 %8.3g S2 %8.3g\n",
                g, nrow(x), if (nrow(u)) max(u$S1) else NA,
                if (nrow(u)) max(u$S2) else NA, nrow(v),
                if (nrow(v)) min(v$S1) else NA, if (nrow(v)) min(v$S2) else NA))
  }
}
cat("\n== listed fits\n")
x <- f[grepl(a[1], f$id), ]
x <- x[order(x$id, x$c), ]
x[c("excess", "S1", "S2")] <- lapply(x[c("excess", "S1", "S2")], signif, 4)
print(x[c("id", "c", "fire", "excess", "S1", "S1par", "S2", "S2par")],
      row.names = FALSE)
