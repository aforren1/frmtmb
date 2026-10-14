# Lane illse: the false-alarm battery's max |SE / ref - 1| (dev/illse-sf-
# rev3-falarm.R), lane against base, per model, on the rows where the
# two differ: how many got closer to the reference, and the largest
# value on each side.
#   Rscript dev/illse-falarm-rel.R <dir> <cfg>
a <- commandArgs(TRUE)
rd <- function(f) utils::read.delim(f, colClasses = "character")
rows <- NULL
for (g in c("nl", "lin", "poly", "mixed", "bounded")) {
  x <- rd(file.path(a[1], sprintf("%s-%s-lane.tsv", g, a[2])))
  y <- rd(file.path(a[1], sprintf("%s-%s-base.tsv", g, a[2])))
  rows <- rbind(rows, data.frame(model = x$model, seed = x$seed,
                                 lane = as.numeric(x$max_rel),
                                 base = as.numeric(y$max_rel)))
}
d <- rows[is.finite(rows$lane) & is.finite(rows$base) &
            rows$lane != rows$base, ]
agg <- do.call(rbind, lapply(split(d, d$model), function(z) {
  data.frame(model = z$model[1], rows = nrow(z),
             closer = sum(z$lane < z$base), lane_max = max(z$lane),
             base_max = max(z$base))
}))
print(agg, row.names = FALSE)
cat("rows that differ:", nrow(d), "| lane closer:", sum(d$lane < d$base), "\n")
