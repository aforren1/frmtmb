# Exact replacement in a tracked source, keeping the checkout's CRLF.
rep1 <- function(path, old, new) {
  s <- readChar(path, file.size(path), useBytes = TRUE)
  if (grepl("\r\n", s, fixed = TRUE)) {
    old <- gsub("\n", "\r\n", old, fixed = TRUE)
    new <- gsub("\n", "\r\n", new, fixed = TRUE)
  }
  hit <- gregexpr(old, s, fixed = TRUE)[[1]]
  if (hit[1] < 0 || length(hit) != 1L) {
    stop("rep1: ", length(hit[hit > 0]), " hits in ", path, " for: ",
         substr(old, 1, 70))
  }
  s <- sub(old, new, s, fixed = TRUE)
  writeBin(charToRaw(s), path)
}
