# Lane optrap: summary of dev/optrap-sweep.sh, per variant and BLAS:
# fits that error, and the seeds where one build errors and the other
# fits.
fs <- Sys.glob("dev/optrap-log/sweep/*.tsv")
X <- lapply(fs, function(f) {
  x <- read.delim(f, quote = "", colClasses = "character")
  p <- strsplit(sub("[.]tsv$", "", basename(f)), "-")[[1]]
  x$variant <- p[1]; x$blas <- p[2]; x$arm <- p[3]
  x
})
X <- do.call(rbind, X)
for (v in unique(X$variant)) for (b in unique(X$blas)) {
  r <- X[X$variant == v & X$blas == b & X$arm == "r7", ]
  l <- X[X$variant == v & X$blas == b & X$arm == "lane", ]
  m <- merge(r, l, by = "seed", suffixes = c(".r7", ".lane"))
  cat(sprintf(paste0("%-5s %-5s seeds %d | errors r7 %d lane %d | r7 errors ",
                     "where lane fits %d | lane errors where r7 fits %d\n"),
              v, b, nrow(m), sum(m$status.r7 == "error"),
              sum(m$status.lane == "error"),
              sum(m$status.r7 == "error" & m$status.lane == "ok"),
              sum(m$status.lane == "error" & m$status.r7 == "ok")))
}
