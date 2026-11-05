# Wrap the review at 80 columns: prose and bullets only, never a table or an indented
# code line. Bullets keep their hanging indent.
f <- "dev/reviews/2026-11-02-predread.md"
x <- readLines(f, warn = FALSE)
out <- character()
for (l in x) {
  if (nchar(l) <= 80 || startsWith(l, "|") || startsWith(l, "    ")) {
    out <- c(out, l)
    next
  }
  lead <- regmatches(l, regexpr("^ *(- )?", l))
  ex <- strrep(" ", nchar(lead))
  body <- substring(l, nchar(lead) + 1L)
  w <- strwrap(body, width = 80 - nchar(lead) + 1L)
  w <- paste0(c(lead, rep(ex, length(w) - 1L)), w)
  out <- c(out, w)
}
writeLines(out, f, useBytes = TRUE)
cat("lines", length(out), " non-table lines over 80:",
    sum(nchar(out) > 80 & !startsWith(out, "|") & !startsWith(out, "    ")), "\n")
