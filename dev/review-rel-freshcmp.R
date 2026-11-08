# Release review 0.80.0: the nlstart reviewer's fresh rows on a candidate
# build against the release build (rellib-r19, dev/rel080-log/rr-final),
# row by row. The best value of a row is the release log's `ref` (best of
# 16-17 started fits) and every default value seen on either build.
#
#   Rscript dev/review-rel-freshcmp.R <candidate dir>
cand <- commandArgs(trailingOnly = TRUE)[1]
rd <- function(files) {
  out <- NULL
  for (f in files) {
    x <- grep("^ROW ", readLines(f, warn = FALSE), value = TRUE)
    p <- strsplit(sub("^ROW ", "", x), " ", fixed = TRUE)
    df <- data.frame(fam = vapply(p, `[`, "", 1), body = vapply(p, `[`, "", 2),
                     seed = vapply(p, `[`, "", 3),
                     ll = suppressWarnings(as.numeric(vapply(p, `[`, "", 4))),
                     nwarn = as.integer(vapply(p, `[`, "", 7)),
                     wclass = vapply(p, `[`, "", 8),
                     multi = vapply(p, `[`, "", 9),
                     ref = suppressWarnings(as.numeric(vapply(p, `[`, "", 13))))
    df$set <- if (grepl("nlfresh2|fresh2-after|/f2[.]txt$", f)) "fresh2" else "fresh"
    out <- rbind(out, df)
  }
  out
}
a2 <- commandArgs(trailingOnly = TRUE)
rel <- rd(if (length(a2) > 1L) Sys.glob(a2[-1L]) else
  Sys.glob("dev/rel080-log/rr-final/nlfresh*/*.txt"))
can <- rd(Sys.glob(file.path(cand, "*.txt")))
k <- c("set", "fam", "body", "seed")
m <- merge(rel, can, by = k, suffixes = c(".rel", ".can"))
cat("rows: release", nrow(rel), " candidate", nrow(can), " matched",
    nrow(m), "\n")
m <- m[is.finite(m$ll.rel) | is.finite(m$ll.can), ]
best <- pmax(m$ref.rel, m$ll.rel, m$ll.can, na.rm = TRUE)
m$short.rel <- best - m$ll.rel
m$short.can <- best - m$ll.can
for (a in c("rel", "can")) {
  s <- m[[paste0("short.", a)]]
  w <- m[[paste0("nwarn.", a)]]
  cat(sprintf("%s: short > 1e-4 %d (silent %d), short > 1 %d (silent %d), multistart %d\n",
              a, sum(s > 1e-4, na.rm = TRUE), sum(s > 1e-4 & w == 0, na.rm = TRUE),
              sum(s > 1, na.rm = TRUE), sum(s > 1 & w == 0, na.rm = TRUE),
              sum(m[[paste0("multi.", a)]] == "1")))
}
d <- m$ll.can - m$ll.rel
cat("candidate - release: higher > 1e-4", sum(d > 1e-4, na.rm = TRUE),
    " lower < -1e-4", sum(d < -1e-4, na.rm = TRUE), "\n")
ch <- m[which(abs(d) > 1e-4 | m$nwarn.rel != m$nwarn.can |
          xor(m$short.rel > 1e-4 & m$nwarn.rel == 0,
              m$short.can > 1e-4 & m$nwarn.can == 0)), ]
print(ch[, c(k, "ll.rel", "nwarn.rel", "multi.rel", "ll.can", "nwarn.can",
             "multi.can", "short.rel", "short.can")], row.names = FALSE)
cat("FRESHCMP DONE\n")
