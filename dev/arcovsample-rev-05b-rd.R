# REVIEW script 05b: the ADDED lines only (the whole-file counts in 05
# include pre-existing text), plus the rendered new sections located by
# a distinctive phrase.
#
#   Rscript dev/arcovsample-rev-05b-rd.R

WT <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample"
outdir <- file.path(WT, "dev/arcovsample-rev-log/rd")

added <- system2("git", c("-C", WT, "diff", "HEAD", "--unified=0"),
                 stdout = TRUE)
added <- grep("^\\+[^+]", added, value = TRUE)
cat("added diff lines: ", length(added), "\n", sep = "")
bad <- c(emdash = "\u2014", endash = "\u2013", lsquo = "\u2018",
         rsquo = "\u2019", ldquo = "\u201c", rdquo = "\u201d")
for (i in seq_along(bad)) {
  k <- grep(bad[i], added, fixed = TRUE, value = TRUE)
  cat("  ", names(bad)[i], ": ", length(k), "\n", sep = "")
  for (s in k) cat("    ", s, "\n", sep = "")
}
sp <- grep(" - ", added, fixed = TRUE, value = TRUE)
cat("  spaced hyphens ' - ': ", length(sp), "\n", sep = "")
for (s in sp) cat("    ", s, "\n", sep = "")
brit <- c("behaviour", "colour", "centre", "organis", "recognis",
          "normalis", "summaris", "analys", "modelling", "labelled",
          "favour", "licence", "practis", "whilst", "learnt", "fulfil",
          "grey")
for (b in brit) {
  k <- grep(b, added, ignore.case = TRUE, value = TRUE)
  if (length(k)) { cat("  BRITISH '", b, "':\n", sep = "")
    for (s in k) cat("    ", s, "\n", sep = "") }
}
long <- added[nchar(added) > 81L]
cat("  added lines over 80 columns (excluding the diff's '+'): ",
    length(long), "\n", sep = "")
for (s in long) cat("    [", nchar(s) - 1L, "] ", s, "\n", sep = "")

cat("\n======== rendered new sections ========\n")
grab <- function(nm, key, nlines = 30L) {
  txt <- readLines(file.path(outdir, nm), warn = FALSE)
  i <- grep(key, txt, fixed = TRUE)
  cat("\n---- ", nm, " at '", key, "' (", length(i), " hit(s))\n",
      sep = "")
  if (!length(i)) return(invisible())
  cat(paste(txt[max(1L, i[1L] - 2L):min(length(txt),
                                        i[1L] + nlines)],
            collapse = "\n"), "\n")
}
grab("sample-log_lik.txt", "Autocorrelation, and what a row conditions")
grab("sample-loo.txt", "A time series, and what is left out", 20L)
grab("frmtmb-sampling-api.txt", "arma_cond_resp(fit)", 20L)
grab("frmtmb-autocor.txt", "s 'log_lik()' and 'loo()' DO cover", 8L)
grab("sample-log_lik.txt", "residual correlation MATRIX", 8L)
cat("\nDONE\n")
