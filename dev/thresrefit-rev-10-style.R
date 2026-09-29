## REVIEW claim 7, style: em dashes, British spellings, spaced hyphens
## standing in for an em dash, 80-column overruns and unescaped percent
## signs, over the files the diff touches only.
files <- c("NEWS.md", "R/bootstrap.R", "R/fit.R", "R/frame.R",
           "R/influence.R", "R/simulate-new.R", "R/thres.R",
           "man/frm.Rd", "man/frm_bootstrap.Rd",
           "man/influence.frmtmb_fit.Rd",
           "tests/testthat/test-thres-refit.R",
           "tests/testthat/test-simulate-ergonomics.R",
           "dev/thresrefit-findings.md")
root <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit"
brit <- paste0("behaviour|colour|centre|centred|centres|analyse|analysed|",
               "organise|organised|normalise|normalised|",
               "modelling|modelled|labelled|labelling|",
               "favour|neighbour|licence|practise|whilst|amongst|",
               "recognise|summarise|utilise|mislabell?ed")
for (f in files) {
  p <- file.path(root, f)
  if (!file.exists(p)) { cat("MISSING", f, "\n"); next }
  x <- readLines(p, warn = FALSE)
  cat("==", f, " lines =", length(x), "\n")
  em <- grep("\u2014|\u2013", x)
  if (length(em)) {
    cat("   EM/EN DASH at lines:", paste(em, collapse = ","), "\n")
    for (i in em) cat("      ", i, ": ", x[i], "\n", sep = "")
  }
  sp <- grep("[[:alnum:]] - [[:alnum:]]", x)
  if (length(sp)) {
    cat("   SPACED HYPHEN at lines:", paste(sp, collapse = ","), "\n")
    for (i in utils::head(sp, 8)) cat("      ", i, ": ", x[i], "\n", sep = "")
  }
  bb <- grep(brit, x, ignore.case = TRUE)
  if (length(bb)) {
    cat("   BRITISH SPELLING candidates at lines:",
        paste(bb, collapse = ","), "\n")
    for (i in bb) cat("      ", i, ": ", x[i], "\n", sep = "")
  }
  emoji <- grep("[\U0001F300-\U0001FAFF\u2600-\u27BF]", x)
  if (length(emoji)) cat("   EMOJI at lines:", paste(emoji, collapse = ","),
                         "\n")
  long <- which(nchar(x) > 80)
  if (length(long)) {
    cat("   OVER 80 COLUMNS at lines:", paste(long, collapse = ","),
        " widths:", paste(nchar(x[long]), collapse = ","), "\n")
  }
  if (grepl("\\.Rd$", f)) {
    pc <- grep("(^|[^\\\\])%", x)
    if (length(pc)) cat("   UNESCAPED % at lines:",
                        paste(pc, collapse = ","), "\n")
  }
}
## the @noRd adjacency rule: a @noRd block must not be followed by text
for (f in c("R/thres.R", "R/frame.R", "R/simulate-new.R", "R/influence.R")) {
  x <- readLines(file.path(root, f), warn = FALSE)
  i <- grep("^#' @noRd\\s*$", x)
  bad <- i[vapply(i, function(k) {
    k < length(x) && grepl("^#'", x[k + 1L]) &&
      !grepl("^#' @", x[k + 1L])
  }, TRUE)]
  cat("== @noRd blocks in", f, ":", length(i),
      " followed by roxygen text:", length(bad),
      if (length(bad)) paste("at", paste(bad, collapse = ",")) else "", "\n")
}
cat("DONE rev-10\n")
