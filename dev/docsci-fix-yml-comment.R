# Rewrap the one comment paragraph this lane edited in the seven
# extension _pkgdown.yml files, so it stays inside 80 columns. Touches
# nothing else: it rewrites only the contiguous run of comment lines that
# holds "build-docs.R", and it refuses a file where that run is not found.
files <- Sys.glob("extensions/*/_pkgdown.yml")
for (f in files) {
  x <- readLines(f, warn = FALSE)
  hit <- grep("build-docs[.]R", x)
  if (length(hit) != 1L) {
    stop("expected one build-docs.R line in ", f, ", found ", length(hit))
  }
  # The paragraph runs from the "# The site builds" line to the first
  # line that is not a comment.
  start <- max(grep("^# The site builds", x))
  end <- start
  while (end < length(x) && startsWith(x[end + 1L], "#")) end <- end + 1L
  body <- paste(trimws(sub("^#[ ]?", "", x[start:end])), collapse = " ")
  wrapped <- paste0("# ", strwrap(body, width = 78))
  writeLines(c(x[seq_len(start - 1L)], wrapped,
               x[seq(end + 1L, length(x))]), f)
  cat(f, ": ", end - start + 1L, " lines -> ", length(wrapped), "\n",
      sep = "")
}
