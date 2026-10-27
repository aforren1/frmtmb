# Fill the generated blocks of a findings draft from the run logs, so
# that no count in the record is typed. Usage:
#   Rscript dev/acprior-fill.R <draft.md> <tag=log ...>
# Each @@TAG@@ line of the draft is replaced by the named log, indented
# four spaces; the result is appended to dev/acprior-findings.md.
a <- commandArgs(TRUE)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
draft <- readLines(a[1])
for (kv in a[-1]) {
  k <- sub("=.*", "", kv)
  v <- sub("^[^=]*=", "", kv)
  blk <- paste0("    ", readLines(v))
  i <- which(draft == paste0("@@", k, "@@"))
  if (length(i) != 1L) stop("no single @@", k, "@@ line")
  draft <- c(draft[seq_len(i - 1L)], blk, draft[-seq_len(i)])
}
if (any(grepl("^@@", draft))) stop("unfilled: ", grep("^@@", draft, value = TRUE))
cat(draft, file = file.path(wt, "dev", "acprior-findings.md"), sep = "\n",
    append = TRUE)
cat("appended", length(draft), "lines\n")
