# Lane adhess: summary of the false-alarm battery's max |SE / ref - 1|,
# lane against base, per model, for one BLAS configuration.
#   Rscript dev/adhess-sf-sum.R <dir with falarm/> <cfg>
a <- commandArgs(TRUE)
rd <- function(f) utils::read.delim(f, colClasses = "character")
out <- list()
for (g in c("nl", "lin", "poly", "mixed", "bounded")) {
  x <- rd(file.path(a[1], "falarm", sprintf("%s-%s-lane.tsv", g, a[2])))
  y <- rd(file.path(a[1], "falarm", sprintf("%s-%s-base.tsv", g, a[2])))
  stopifnot(identical(x$model, y$model), identical(x$seed, y$seed))
  out[[g]] <- data.frame(group = g, model = x$model, seed = x$seed,
                         lane = suppressWarnings(as.numeric(x$max_rel)),
                         base = suppressWarnings(as.numeric(y$max_rel)))
}
d <- do.call(rbind, out)
ch <- !(is.na(d$lane) & is.na(d$base)) &
  (is.na(d$lane) != is.na(d$base) | d$lane != d$base)
cat(sprintf("rows %d | max_rel changed %d | lane lower %d | lane higher %d | NA on one arm %d\n",
            nrow(d), sum(ch, na.rm = TRUE),
            sum(ch & d$lane < d$base, na.rm = TRUE),
            sum(ch & d$lane > d$base, na.rm = TRUE),
            sum(is.na(d$lane) != is.na(d$base))))
s <- do.call(rbind, lapply(split(d, d$model), function(m) {
  data.frame(model = m$model[1], group = m$group[1], n = nrow(m),
             changed = sum(m$lane != m$base, na.rm = TRUE),
             base_max = max(m$base, na.rm = TRUE),
             lane_max = max(m$lane, na.rm = TRUE),
             lane_higher = sum(m$lane > m$base, na.rm = TRUE))
}))
s <- s[order(s$group, s$model), ]
print(format(s, digits = 3), row.names = FALSE)
hi <- d[ch & d$lane > d$base & !is.na(ch), ]
if (nrow(hi)) {
  cat("\nrows where the lane is farther from the reference:\n")
  print(format(hi, digits = 4), row.names = FALSE)
}
