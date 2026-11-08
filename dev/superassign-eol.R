# Report, per file, how many lines end in CRLF and how many in a bare
# LF, and every line over 80 characters. With `--fix`, rewrite a file
# that mixes the two so every line ends in CRLF.
#
# Usage: Rscript dev/superassign-eol.R [--fix] <file>...

a <- commandArgs(trailingOnly = TRUE)
fix <- "--fix" %in% a
for (f in setdiff(a, "--fix")) {
  txt <- rawToChar(readBin(f, "raw", file.info(f)$size))
  crlf <- lengths(regmatches(txt, gregexpr("\r\n", txt, fixed = TRUE)))
  lf <- lengths(regmatches(txt, gregexpr("\n", txt, fixed = TRUE))) - crlf
  lines <- strsplit(gsub("\r\n", "\n", txt, fixed = TRUE), "\n",
                    fixed = TRUE)[[1L]]
  long <- which(nchar(lines) > 80L)
  cat(sprintf("%s: crlf=%d lf=%d long=%s\n", f, crlf, lf,
              if (length(long)) paste(long, collapse = ",") else "-"))
  if (fix && crlf > 0L && lf > 0L) {
    out <- paste0(paste(lines, collapse = "\r\n"), "\r\n")
    writeBin(charToRaw(out), f)
    cat("  rewritten as CRLF\n")
  }
}
