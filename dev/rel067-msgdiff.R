# Rows of the ported brms tier whose recorded message moved between the
# 0.66.0 ledger and this round's record, with their verdict and reason.
# A moved message on a row that still does not hold means its recorded
# reason may be stale; each is read by hand at consolidation.
#
#   Rscript dev/rel067-msgdiff.R > dev/rel067-log/msgdiff.txt
rd <- function(f, ...) utils::read.delim(f, quote = "",
                                         colClasses = "character",
                                         na.strings = NULL, ...)
old <- rd("dev/rel067-ledger-before.tsv")
recs <- list.files("dev/brmsport-log", "^rec-.*[.]tsv$", full.names = TRUE)
rec <- do.call(rbind, lapply(recs, function(f) {
  x <- rd(f, header = FALSE)
  data.frame(pkg = x[[2]], id = x[[3]], verdict = x[[4]], msg = x[[7]])
}))
man <- rbind(rd("dev/brmsport-verdicts-manual.tsv"),
             rd("dev/brmsport-verdicts-manual-fit.tsv"))
norm <- function(s) gsub("[[:space:]]+", " ", gsub("0x[0-9a-f]+", "", s))
n <- 0L
for (i in seq_len(nrow(rec))) {
  r <- rec[i, ]
  o <- old[old$id == r$id, , drop = FALSE]
  if (!nrow(o) || all(o$outcome == "pass")) next
  om <- norm(o$message[1]); nm <- norm(r$msg)
  if (identical(om, nm)) next
  m <- man[man$id == r$id, , drop = FALSE]
  n <- n + 1L
  cat("==", r$id, r$pkg, "| ledger:", o$outcome[1], "/", o$class[1], "\n")
  cat("   reason: ", if (nrow(m)) substr(m$reason[1], 1, 300) else "-", "\n")
  cat("   before: ", substr(om, 1, 300), "\n")
  cat("   now:    ", substr(nm, 1, 300), "\n")
}
cat("rows with a moved message:", n, "\n")
