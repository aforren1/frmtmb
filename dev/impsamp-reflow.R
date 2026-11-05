# Rewrap the paragraphs of the findings that pass 80 columns; tables, code
# and list markers are left as they are.
p <- "dev/impsamp-findings.md"
x <- readLines(p)
out <- character(0)
i <- 1L
while (i <= length(x)) {
  ln <- x[[i]]
  if (!nzchar(ln) || grepl("^(#|\\||    )", ln)) {
    out <- c(out, ln)
    i <- i + 1L
    next
  }
  j <- i
  while (j < length(x) && nzchar(x[[j + 1L]]) &&
           !grepl("^(#|\\||    |- )", x[[j + 1L]])) {
    j <- j + 1L
  }
  blk <- x[i:j]
  if (any(nchar(blk) > 80L)) {
    bullet <- grepl("^- ", blk[[1L]])
    blk <- strwrap(paste(blk, collapse = " "), width = 79,
                   exdent = if (bullet) 2L else 0L)
  }
  out <- c(out, blk)
  i <- j + 1L
}
writeLines(out, p)
cat("lines over 80, outside tables:",
    sum(nchar(out) > 80L & !grepl("^\\|", out)), "\n")
