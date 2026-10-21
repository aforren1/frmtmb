# The re-measure of rellib-r10 (tag r10) repeats the 0.72.0
# consolidation's r072 run of the same release. Compare r10's columns
# of each outcome table with r072's; where every one agrees, drop r10's
# columns and its row of dev/brms-port/outcome-builds.tsv, so that each
# build has one column. Stops on a difference.
#   Rscript dev/rel073-outcomes-r10.R
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL, check.names = FALSE)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
for (f in c("dev/brms-port/outcomes-hand.tsv",
            "dev/brms-port/outcomes-mech.tsv")) {
  x <- rd(f)
  r10 <- grep("^r10[.]", names(x), value = TRUE)
  if (!length(r10)) { cat(f, ": no r10 column\n"); next }
  for (k in r10) {
    m <- sub("^r10[.]", "r072.", k)
    stopifnot(m %in% names(x))
    cat(f, ":", k, "==", m, "on", sum(x[[k]] == x[[m]]), "of", nrow(x),
        "rows\n")
    stopifnot(identical(x[[k]], x[[m]]))
  }
  x[r10] <- NULL
  wr(x, f)
  cat("  dropped", r10, "\n")
}
b <- readLines("dev/brms-port/outcome-builds.tsv")
b <- b[!startsWith(b, "r10\t")]
writeLines(b, "dev/brms-port/outcome-builds.tsv")
