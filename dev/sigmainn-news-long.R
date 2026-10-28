# Sentences over 25 words in the development sections of the NEWS files
# this lane writes (a code span counts as one word), lines over 80, em
# dashes. Usage: Rscript dev/sigmainn-news-long.R
for (f in c("NEWS.md", "extensions/frmtmb.sample/NEWS.md")) {
  x <- readLines(f)
  end <- grep("^# ", x)[2L] - 1L
  sec <- x[seq_len(end)]
  cat(f, ": lines over 80:", sum(nchar(sec) > 80),
      "; em dashes:", sum(grepl("—", sec)), "\n")
  txt <- paste(sub("^[*] ", "", trimws(sec[!grepl("^#", sec)])),
               collapse = " ")
  txt <- gsub("`[^`]*`", "CODE", txt)
  s <- unlist(strsplit(txt, "(?<=[.:])\\s+(?=[A-Z*(])", perl = TRUE))
  n <- lengths(strsplit(trimws(s), "\\s+"))
  long <- s[n > 25]
  cat("  sentences over 25 words:", length(long), "\n")
  for (l in long) cat("   -", substr(l, 1, 160), "\n")
}
