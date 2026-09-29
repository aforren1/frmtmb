# REVIEW script 05, claim 6: render the changed Rd topics and grep the
# RENDERED text, not the source. `%` starts a comment in Rd even inside
# a verbatim macro, so a source read cannot settle this.
#
#   Rscript dev/arcovsample-rev-05-rd.R

WT <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample"
outdir <- file.path(WT, "dev/arcovsample-rev-log/rd")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  file.path(WT, "extensions/frmtmb.sample/man/sample-log_lik.Rd"),
  file.path(WT, "extensions/frmtmb.sample/man/sample-loo.Rd"),
  file.path(WT, "man/frmtmb-sampling-api.Rd"),
  file.path(WT, "man/frmtmb-autocor.Rd"))

# U+2014 em dash, U+2013 en dash, U+2018/9/C/D curly quotes: named by
# code point so that this file's own encoding cannot hide one.
bad_chars <- c(emdash = "\u2014", endash = "\u2013",
               lsquo = "\u2018", rsquo = "\u2019",
               ldquo = "\u201c", rdquo = "\u201d")
# British spellings that the house style bans. Word boundaries, so
# "colourless" and "behaviourism" are caught too.
brit <- c("behaviour", "colour", "centre", "organise", "organisation",
          "recognise", "normalise", "summarise", "analyse", "modelling",
          "labelled", "favour", "licence", "practise", "whilst",
          "towards", "learnt", "fulfil", "grey")

for (f in files) {
  nm <- basename(f)
  txt <- capture.output(tools::Rd2txt(f, out = stdout(), package = "x"))
  writeLines(txt, file.path(outdir, sub("\\.Rd$", ".txt", nm)))
  one <- paste(txt, collapse = "\n")
  cat("\n================ ", nm, " (", length(txt),
      " rendered lines)\n", sep = "")
  for (i in seq_along(bad_chars)) {
    k <- grep(bad_chars[i], txt, fixed = TRUE)
    cat("  ", names(bad_chars)[i], ": ", length(k), "\n", sep = "")
  }
  hits <- character(0)
  for (b in brit) {
    if (grepl(b, one, ignore.case = TRUE)) hits <- c(hits, b)
  }
  cat("  British spellings: ",
      if (length(hits)) paste(hits, collapse = ", ") else "none", "\n",
      sep = "")
  # " - " standing in for an em dash
  cat("  ' - ' spaced hyphens: ", sum(grepl(" - ", txt, fixed = TRUE)),
      "\n", sep = "")
  cat("  lines over 80 columns: ", sum(nchar(txt) > 80), "\n", sep = "")
}

cat("\n\n======== the new sections, as rendered ========\n")
show <- function(nm, from, to) {
  txt <- readLines(file.path(outdir, nm))
  i <- grep(from, txt)[1L]
  j <- if (is.na(to)) length(txt) else grep(to, txt)[1L]
  cat("\n---- ", nm, ": ", from, "\n", sep = "")
  cat(paste(txt[i:max(i, j - 1L)], collapse = "\n"), "\n")
}
show("sample-log_lik.txt", "Autocorrelation, and what a row conditions on",
     "Multivariate models")
show("sample-loo.txt", "A time series, and what is left out",
     "Priors, and what these numbers mean")

cat("\n\n======== NEWS ========\n")
for (f in c(file.path(WT, "NEWS.md"),
            file.path(WT, "extensions/frmtmb.sample/NEWS.md"))) {
  txt <- readLines(f, warn = FALSE)
  h <- grep("^# ", txt)
  cat("\n---- ", f, "\n", sep = "")
  cat("first heading: ", txt[h[1L]], "\n", sep = "")
  cat("second heading: ", txt[h[2L]], "\n", sep = "")
  blk <- txt[h[1L]:(h[2L] - 1L)]
  cat("  em dashes: ", sum(grepl("\u2014", blk, fixed = TRUE)),
      "  ' - ': ", sum(grepl(" - ", blk, fixed = TRUE)),
      "  over 80 cols: ", sum(nchar(blk) > 80), "\n", sep = "")
  # a version number inside the development-version block
  vn <- grep("[0-9]+[.][0-9]+[.][0-9]+", blk, value = TRUE)
  cat("  lines naming a version number: ", length(vn), "\n", sep = "")
  for (v in vn) cat("    ", v, "\n", sep = "")
  cat("  --- block ---\n")
  cat(paste(blk, collapse = "\n"), "\n")
}
cat("\nDONE\n")
