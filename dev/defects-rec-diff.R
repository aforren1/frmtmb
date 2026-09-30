# Lane wt-defects: rows whose recorded outcome or message differs between
# two sets of rec-*.tsv files.
#   Rscript dev/defects-rec-diff.R <dir A> <dir B>
a <- commandArgs(trailingOnly = TRUE)
cols <- c("kind", "pkg", "id", "verdict", "held", "vacuous", "msg",
          "caught", "raw_held")
read_dir <- function(d) {
  fs <- list.files(d, pattern = "^rec-.*[.]tsv$", full.names = TRUE)
  x <- do.call(rbind, lapply(fs, function(f) {
    utils::read.delim(f, header = FALSE, quote = "", col.names = cols,
                      colClasses = "character", na.strings = NULL)
  }))
  x <- x[x$kind == "assert", ]
  x$key <- paste(x$pkg, x$id)
  x
}
A <- read_dir(a[1])
B <- read_dir(a[2])
m <- merge(A[c("key", "held", "msg")], B[c("key", "held", "msg")],
           by = "key", suffixes = c(".a", ".b"), all = TRUE)
# a number drawn at random differs run to run, so compare messages with
# the digits removed
strip <- function(s) gsub("[0-9.e+-]+", "#", s)
chg <- m[is.na(m$held.a) | is.na(m$held.b) | m$held.a != m$held.b |
           strip(m$msg.a) != strip(m$msg.b), ]
for (i in seq_len(nrow(chg))) {
  cat(sprintf("%s  held %s -> %s\n   A: %s\n   B: %s\n", chg$key[i],
              chg$held.a[i], chg$held.b[i], substr(chg$msg.a[i], 1, 200),
              substr(chg$msg.b[i], 1, 200)))
}
cat("rows compared:", nrow(m), " changed:", nrow(chg), "\n")
