# Append the section to the review, wrapped at 80 columns.
d <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev"
s <- readLines(file.path(d, "review-impsamp-fr2-section.txt"))
out <- character(0)
k <- 1L
n <- length(s)
while (k <= n) {
  ln <- s[k]
  if (!nzchar(ln) || grepl("^(\\||#)", ln)) {
    out <- c(out, ln)
    k <- k + 1L
    next
  }
  e <- k
  while (e < n && nzchar(s[e + 1L]) && !grepl("^(\\||#|- )", s[e + 1L])) e <- e + 1L
  txt <- paste(trimws(s[k:e]), collapse = " ")
  b <- grepl("^- ", ln)
  out <- c(out, strwrap(sub("^- ", "", txt), width = 79, initial = if (b) "- " else "",
                        prefix = if (b) "  " else ""))
  k <- e + 1L
}
p <- file.path(d, "reviews/2026-11-02-impsamp.md")
stopifnot(!any(grepl("Short check of fix round 2", readLines(p), fixed = TRUE)))
cat(out, file = p, sep = "\n", append = TRUE)
cat("appended", length(out), "lines\n")
