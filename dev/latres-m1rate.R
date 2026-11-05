# Hit and false-alarm rates of the M1 warning over dev/latres-m1.R's data sets.
# A design is bad when the median sderr ratio to the exact MLE leaves [0.9, 1.15].
# Rule: median information below `floor`, or at most `tr` trials.
fs <- list.files("dev/latres-out/m1", "[.]tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.delim(f, stringsAsFactors = FALSE)))
x$ratio <- x$sderr / pmin(x$x_sderr, 10)
tr <- c(binom3 = 3)
ds <- split(x, x$design)
sumd <- do.call(rbind, lapply(ds, function(z) data.frame(
  design = z$design[1], n = nrow(z), ratio = median(z$ratio, na.rm = TRUE),
  ar = median(z$ar, na.rm = TRUE), x_ar = median(z$x_ar, na.rm = TRUE))))
sumd$bad <- sumd$ratio < 0.9 | sumd$ratio > 1.15 | abs(sumd$ar - sumd$x_ar) > 0.15
for (fl in c(0.3, 0.35, 0.4, 0.45, 0.5, 0.6)) {
  warn <- x$info_med < fl | x$design %in% names(tr)
  w <- tapply(warn, x$design, mean)
  sumd[[paste0("w", fl)]] <- w[sumd$design]
  cat(sprintf("floor %.2f: hit %.3f (%d data sets), false alarm %.3f (%d)\n", fl,
              mean(warn[x$design %in% sumd$design[sumd$bad]]),
              sum(x$design %in% sumd$design[sumd$bad]),
              mean(warn[x$design %in% sumd$design[!sumd$bad]]),
              sum(x$design %in% sumd$design[!sumd$bad])))
}
print(format(sumd[order(sumd$ratio), ], digits = 3), row.names = FALSE)
