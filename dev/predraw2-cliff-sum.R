# Summary of dev/predraw2-cliff-cov.R: coverage by map and by the narrowed
# side's ratio r, per design, and over all fits (the fits with no narrowed
# side draw the same under every map).
# Usage: Rscript dev/predraw2-cliff-sum.R <tsv>...
d <- do.call(rbind, lapply(commandArgs(TRUE), utils::read.delim))
d$design <- sub(":.*", "", d$tag)
arms <- c("cur", "normal", "blend0", "blend25", "blend40")
for (g in unique(d$design)) {
  x <- d[d$design == g, ]
  ch <- x[!is.na(x$r), ]
  cat(sprintf("\n%s: %d fits, %d with a narrowed side\n", g,
              length(unique(x$tag)), length(unique(ch$tag))))
  bins <- cut(ch$r, c(0, 0.1, 0.25, 0.4, 0.5, 1), include.lowest = TRUE)
  for (b in levels(bins)) {
    y <- ch[bins == b, ]
    if (!nrow(y)) next
    cat(sprintf("  r %-10s n %3d", b, length(unique(y$tag))))
    for (a in arms) {
      z <- y[y$arm == a, ]
      cat(sprintf(" | %s %.3f/%.3f", a, mean(z$cov90), mean(z$cov98)))
    }
    cat("\n")
  }
  for (a in arms) {
    z <- rbind(x[is.na(x$r), ], x[!is.na(x$r) & x$arm == a, ])
    cat(sprintf("  all fits, %-8s cov90 %.3f cov98 %.3f median w98 %.3f\n", a,
                mean(z$cov90), mean(z$cov98), stats::median(z$w98)))
  }
}
# pooled over the designs with a known law (R and S), by bin of r
ch <- d[!is.na(d$r) & d$design != "cliff", ]
if (nrow(ch)) {
  cat("\npooled over R and S designs, fits with a narrowed side:\n")
  bins <- cut(ch$r, c(0, 0.1, 0.25, 0.4, 0.5, 1), include.lowest = TRUE)
  for (b in levels(bins)) {
    y <- ch[bins == b, ]
    if (!nrow(y)) next
    cat(sprintf("  r %-10s n %3d", b, length(unique(y$tag))))
    for (a in arms) {
      z <- y[y$arm == a, ]
      if (nrow(z)) {
        cat(sprintf(" | %s %.4f/%.4f", a, mean(z$cov90), mean(z$cov98)))
      }
    }
    cat("\n")
  }
}
