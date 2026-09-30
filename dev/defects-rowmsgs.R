# Lane wt-defects: the recorded outcome of every row the 0.65.0 ledger
# (dev/defects-ledger-before.tsv) calls a defect, read from the current
# dev/brmsport-log/rec-*.tsv files.
#   Rscript dev/defects-rowmsgs.R [id ...]
recs <- list.files("dev/brmsport-log", pattern = "^rec-.*[.]tsv$",
                   full.names = TRUE)
cols <- c("kind", "pkg", "id", "verdict", "held", "vacuous", "msg",
          "caught", "raw_held")
rec <- do.call(rbind, lapply(recs, function(f) {
  x <- utils::read.delim(f, header = FALSE, quote = "", col.names = cols,
                         colClasses = "character", na.strings = NULL)
  x[x$kind == "assert", ]
}))
before <- utils::read.delim("dev/defects-ledger-before.tsv", quote = "",
                            colClasses = "character", na.strings = NULL)
ids <- commandArgs(trailingOnly = TRUE)
if (!length(ids)) ids <- before$id[before$outcome == "defect"]
for (id in ids) {
  r <- rec[rec$id == id, ]
  for (k in seq_len(nrow(r))) {
    cat(sprintf("%-26s %-13s held=%s  %s\n", id, r$pkg[k], r$held[k],
                substr(r$msg[k], 1, 230)))
  }
}
