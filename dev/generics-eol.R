# Put back the LF line endings my edit scripts destroyed, on the files
# THIS LANE OWNS and on nothing else.
#
# What went wrong, twice, and it is worth writing down.
#
# 1. `writeLines(x, path)` opens a TEXT connection, and on Windows a
#    text connection translates "\n" to "\r\n".  Every file
#    dev/generics-edit*.R rewrote came out CRLF in a tree that is LF,
#    which git reports as "CRLF will be replaced by LF the next time
#    Git touches it".  Found by an edit failing to match its own
#    pattern because of the stray "\r".
#
# 2. The FIRST repair was worse than the defect.  It walked every
#    tracked file, and its condition, "contains a CR followed by an
#    LF", is true of BINARY files too: a gzip byte of 0x0D beside 0x0A
#    is ordinary.  It rewrote 123 `.rds`, `.png`, `.ttf` and `.woff2`
#    files with those bytes removed, corrupting all of them, and said
#    "converted 424 tracked files back to LF" as though that were a
#    result.  They were restored from the index with
#    `git restore --worktree` and re-read to confirm, and this version
#    is the narrow one: an explicit list of the lane's own files, and a
#    refusal on anything that is not plain text.
#
# This is the project's own rule about guards, in a place nobody
# thought to look: the guard's condition must be one that is true ONLY
# of the thing being guarded, and "construct the case where the
# guarded thing is absent" means running it on a file it must not
# touch.
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
own <- c("DESCRIPTION", "NAMESPACE", "_pkgdown.yml", "NEWS.md",
         "R/draws-generics.R", "R/generic-owners.R", "R/loo.R",
         "R/methods-fit.R", "R/predict.R", "R/scales.R", "R/sugar.R",
         "R/zzz.R",
         "tests/testthat/test-generic-collision.R",
         "tests/testthat/test-scale-contract.R",
         "extensions/frmtmb.sample/R/methods-draws.R",
         "dev/generics-findings.md")
own <- c(own, list.files(file.path(root, "dev"), pattern = "^generics-",
                         recursive = TRUE, full.names = FALSE,
                         all.files = FALSE))
own <- unique(own)

is_text <- function(raw) {
  if (!length(raw)) return(TRUE)
  # a NUL byte is the cheap, reliable signature of a binary file, and
  # every source file this lane writes is ASCII
  !any(raw == as.raw(0L)) && all(raw < as.raw(128L))
}

fixed <- character()
skipped <- character()
for (f in own) {
  p <- file.path(root, if (grepl("^generics-", f)) file.path("dev", f)
                 else f)
  if (!file.exists(p) || dir.exists(p)) next
  raw <- readBin(p, "raw", file.info(p)$size)
  if (!is_text(raw)) { skipped <- c(skipped, f); next }
  crlf <- which(raw == as.raw(13L))
  crlf <- crlf[crlf < length(raw) & raw[pmin(crlf + 1L, length(raw))] ==
                 as.raw(10L)]
  if (!length(crlf)) next
  writeBin(raw[-crlf], p)
  fixed <- c(fixed, f)
}
cat(sprintf("converted %d of the lane's own files to LF\n", length(fixed)))
for (f in fixed) cat("  ", f, "\n")
if (length(skipped)) {
  cat("refused (not plain text):\n")
  for (f in skipped) cat("  ", f, "\n")
}
