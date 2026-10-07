# Does seeding the ported brms blocks change any assertion's outcome?
# Compares the "held" column of the base (unseeded) and lane (seeded)
# records written by dev/ciharden-brmsstable.sh.
rd <- function(d) {
  fs <- list.files(d, "^rec-.*[.]tsv$", full.names = TRUE)
  x <- do.call(rbind, lapply(fs, function(f) {
    utils::read.delim(f, header = FALSE, quote = "", colClasses = "character",
                      na.strings = NULL)
  }))
  data.frame(key = paste(x[[1]], x[[2]], x[[3]]), verdict = x[[4]],
             held = x[[5]])
}
b <- rd("dev/ciharden-log/brmsstable-base/A")
l <- rd("dev/ciharden-log/brmsstable-lane/A")
m <- merge(b, l, by = "key", suffixes = c(".base", ".lane"))
cat("rows:", nrow(b), nrow(l), "matched", nrow(m), "\n")
d <- m[m$held.base != m$held.lane, ]
cat("outcome changed by seeding:", nrow(d), "\n")
print(d, row.names = FALSE)
