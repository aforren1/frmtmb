# Append the section to the review, wrapped at 80 columns, with the mutant's count
# read from its log.
d <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev"
s <- readLines(file.path(d, "review-impsamp-fr1-section.txt"))
lg <- readLines(file.path(d, "review-impsamp-out/fr1/mut/test-fix_tau.log"))
r <- grep("^RESULT", lg, value = TRUE)
stopifnot(length(r) == 1L)
s <- sub("FIXTAU_RESULT", sub("^RESULT fix_tau ", "test-importance.R on it: ", r), s, fixed = TRUE)
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
stopifnot(!any(grepl("Short check of fix round 1", readLines(p), fixed = TRUE)))
cat(out, file = p, sep = "\n", append = TRUE)
cat("appended", length(out), "lines; result:", r, "\n")
