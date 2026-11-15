# Rewrap lane quadseam2's CLOSED lines in dev/test-backlog.md to 80
# columns (dev/quadseam2-backlog-edit.R wrote the first line long).
f <- "dev/test-backlog.md"
x <- readLines(f)
at <- grep("^  CLOSED in the round of 2026-11-13 [(]lane quadseam2[)]", x)
for (i in rev(at)) {
  j <- i
  while (j + 1L <= length(x) && grepl("^  [^ -]", x[j + 1L]) &&
         !grepl("^  CLOSED", x[j + 1L])) j <- j + 1L
  txt <- paste(trimws(x[i:j]), collapse = " ")
  w <- strwrap(txt, width = 78, indent = 2, exdent = 2)
  x <- c(x[seq_len(i - 1L)], w, x[-seq_len(j)])
}
writeLines(x, f)
