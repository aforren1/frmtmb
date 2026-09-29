## Verify the touched Rd by RENDERING it, not by reading the source, and
## scan the rendered text for the style rules a reviewer checks.
.libPaths(c("C:/Users/adf44/source/r/wt-thresrefit-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit"
marks <- c("Each deletion refits", "A refit inside the package",
           "Every refit reuses", "A refit that fails",
           "Two deletions cannot", "is not a refit",
           "Calling .frm.. again")
files <- c("influence.frmtmb_fit.Rd", "frm.Rd", "frm_bootstrap.Rd")
bad <- character(0)
for (f in files) {
  out <- tempfile()
  tools::Rd2txt(file.path(root, "man", f), out = out)
  txt <- readLines(out, warn = FALSE)
  cat("=====", f, "=====\n")
  hit <- sort(unique(unlist(lapply(marks, function(m) grep(m, txt)))))
  for (k in hit) {
    # print the whole paragraph the mark sits in
    a <- k
    while (a > 1L && nzchar(trimws(txt[a - 1L]))) a <- a - 1L
    b <- k
    while (b < length(txt) && nzchar(trimws(txt[b + 1L]))) b <- b + 1L
    cat(paste(txt[a:b], collapse = "\n"), "\n\n")
  }
  # the style scan runs on the RENDERED text, which is what a reader sees
  raw <- paste(txt, collapse = "\n")
  for (pat in c("\u2014", "\u2013", "neighbour", "behaviour", "colour",
                "modelling", "mislabelled", "organis", "analyse")) {
    if (grepl(pat, raw, fixed = TRUE)) bad <- c(bad, paste0(f, ": ", pat))
  }
  # Rd2txt underlines a heading the way a teletype did, with an
  # UNDERSCORE, A BACKSPACE and the character: "_\bT_\bh_\be". Both the
  # underscore and the backspace have to go, or the width being measured
  # is still markup. Stripping only the underscore left the backspaces
  # behind, which reported the 46-column heading "The Laplace
  # approximation, and how to check it:" as 85 columns and made a
  # `grepl()` for its own text fail. That produced one phantom
  # over-80 line in this lane's earlier runs, and an exemption written to
  # excuse it.
  plain <- gsub("_\b", "", txt, useBytes = TRUE)
  wide <- which(nchar(plain) > 80)
  if (length(wide)) {
    bad <- c(bad, paste0(f, ": ", length(wide), " rendered lines over 80 (",
                         paste(wide, collapse = ","), "): ",
                         paste0("\"", trimws(plain[wide]), "\"",
                                collapse = ", ")))
  }
}
cat("style problems:",
    if (length(bad)) paste(bad, collapse = "; ") else "none", "\n")
