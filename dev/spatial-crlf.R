# Write the named files with CRLF line endings, as a checkout with
# core.autocrlf = true writes them: every LF becomes CRLF, and a last
# line without a newline is left without one (sed 's/$/\r/' would give it
# a lone CR, which git then reads as binary). Files already CRLF are left
# alone. Usage: Rscript dev/spatial-crlf.R <file> ...
for (f in commandArgs(trailingOnly = TRUE)) {
  x <- readBin(f, "raw", file.size(f))
  if (any(x == as.raw(13))) { cat("skip (has CR):", f, "\n"); next }
  lf <- which(x == as.raw(10))
  out <- raw(length(x) + length(lf))
  pos <- seq_along(x) + findInterval(seq_along(x), lf)
  out[pos] <- x
  out[lf + seq_along(lf) - 1L] <- as.raw(13)
  writeBin(out, f)
  cat("crlf:", f, "\n")
}
