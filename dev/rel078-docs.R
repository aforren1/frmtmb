# Section 6 of dev/rel078-resolve.R alone: dev/test-backlog.md and
# dev/round-handoff.md from their HEAD blobs with this round's sections
# (dev/rel078-backlog-section.md, dev/rel078-handoff-section.md)
# inserted, the same text a full merge rerun writes, stored CRLF. For an
# edit of the two sections after the last merge run, which would
# otherwise need a rerun of the whole merge and a reinstall.
#   Rscript dev/rel078-docs.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel078")
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
head_blob <- function(f) {
  sub("\r$", "", system2("git", c("show", paste0("HEAD:", f)), stdout = TRUE))
}
ins <- function(x, old, new) {
  s <- paste(x, collapse = "\n")
  k <- sum(gregexpr(old, s, fixed = TRUE)[[1]] > 0)
  if (k != 1L) stop("expected 1 match, found ", k, ": ", substr(old, 1, 60))
  strsplit(sub(old, new, s, fixed = TRUE), "\n", fixed = TRUE)[[1]]
}
wr <- function(x, f) {
  con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
  cat("wrote", f, "\n")
}
bl <- ins(head_blob("dev/test-backlog.md"), "\n## Reference\n",
          paste0("\n", paste(rd("dev/rel078-backlog-section.md"),
                             collapse = "\n"), "\n\n## Reference\n"))
wr(bl, "dev/test-backlog.md")
hf <- ins(head_blob("dev/round-handoff.md"),
"0.76.0 release and 2026-10-28 at the 0.77.0 release. Read this, then",
"0.76.0 release, 2026-10-28 at the 0.77.0 release and 2026-10-29 at the
0.78.0 release. Read this, then")
hf <- ins(hf, "\n## The 0.77.0 round, ",
          paste0("\n", paste(rd("dev/rel078-handoff-section.md"),
                             collapse = "\n"), "\n\n## The 0.77.0 round, "))
wr(hf, "dev/round-handoff.md")
