# The re-measure of rellib-r9 (tag r9) repeats the 0.71.0 consolidation's
# "merged" run of the same library. Compare r9's columns of each outcome
# table with merged's; where every one agrees, drop r9's columns and its
# row of dev/brms-port/outcome-builds.tsv, so that each build has one
# column. Stops on a difference.
#   Rscript dev/rel072-outcomes-r9.R
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL, check.names = FALSE)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
for (f in c("dev/brms-port/outcomes-hand.tsv",
            "dev/brms-port/outcomes-mech.tsv")) {
  x <- rd(f)
  r9 <- grep("^r9[.]", names(x), value = TRUE)
  if (!length(r9)) { cat(f, ": no r9 column\n"); next }
  for (k in r9) {
    m <- sub("^r9[.]", "merged.", k)
    stopifnot(m %in% names(x))
    cat(f, ":", k, "==", m, "on", sum(x[[k]] == x[[m]]), "of", nrow(x),
        "rows\n")
    stopifnot(identical(x[[k]], x[[m]]))
  }
  x[r9] <- NULL
  wr(x, f)
  cat("  dropped", r9, "\n")
}
b <- readLines("dev/brms-port/outcome-builds.tsv")
b <- b[!startsWith(b, "r9\t")]
writeLines(b, "dev/brms-port/outcome-builds.tsv")
