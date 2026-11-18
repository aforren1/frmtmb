# Lane drawdir: the largest joint excesses (S2 of dev/drawdir-an2.R) at
# c = 2 over the probe2 files given, with the parameter and the side.
# Usage: Rscript dev/drawdir-top.R <n> <regex of ids> <tsv> [...]
a <- commandArgs(TRUE)
d <- do.call(rbind, lapply(a[-(1:2)], utils::read.delim, stringsAsFactors = FALSE))
d <- d[d$c == 2 & grepl(a[2], d$id), ]
d$S2p <- d$eAp - pmax(d$eBp, 0)
d$S2m <- d$eAm - pmax(d$eBm, 0)
d$S2 <- pmax(d$S2p, d$S2m)
f <- do.call(rbind, lapply(split(d, d$id), function(x) {
  j <- which.max(x$S2)
  data.frame(id = x$id[1], fire = x$fire[1], np = x$np[1], S2 = x$S2[j],
             par = x$par[j], side = if (x$S2p[j] >= x$S2m[j]) "+" else "-",
             eA = max(x$eAp[j], x$eAm[j]), sratio = x$ratio[j])
}))
f <- f[order(-f$S2), ]
cat("fits", nrow(f), "| S2 > 6:", sum(f$S2 > 6, na.rm = TRUE), "| > 10.5:",
    sum(f$S2 > 10.5, na.rm = TRUE), "| NA:", sum(is.na(f$S2)), "\n")
f[c("S2", "eA", "sratio")] <- lapply(f[c("S2", "eA", "sratio")], signif, 4)
print(utils::head(f, as.integer(a[1])), row.names = FALSE)
