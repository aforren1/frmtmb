# Emit the R CMD check status lines for both checked packages as a
# block.  Read from 00check.log, not from the console transcript, so a
# truncated transcript cannot produce a cleaner-looking result than the
# check recorded.
paths <- c(
  frmtmb = "C:/Users/adf44/source/r/generics-check-core/frmtmb.Rcheck/00check.log",
  frmtmb.sample = paste0("C:/Users/adf44/source/r/generics-check-sample/",
                         "frmtmb.sample.Rcheck/00check.log")
)
cat("```\n")
cat("== R CMD check --as-cran, built WITH vignettes, no --no-manual ==\n")
cat("dev/generics-cran.sh, mirroring dev/release/run-check.ps1\n\n")
for (p in names(paths)) {
  f <- paths[[p]]
  if (!file.exists(f)) { cat(sprintf("%-14s NO LOG\n", p)); next }
  ln <- readLines(f, warn = FALSE)
  st <- grep("^Status:", ln, value = TRUE)
  cat(sprintf("%-14s %s\n", p, if (length(st)) st[1] else "NO STATUS"))
  hits <- grep("^\\* checking .*(NOTE|WARNING|ERROR)$", ln)
  for (i in hits) {
    cat(sprintf("   %s\n", sub("^\\* checking ", "", ln[i])))
    j <- i + 1L
    while (j <= length(ln) && !grepl("^\\* ", ln[j]) && nzchar(ln[j])) {
      cat(sprintf("     %s\n", ln[j]))
      j <- j + 1L
    }
  }
  cat("\n")
}
cat("```\n")
