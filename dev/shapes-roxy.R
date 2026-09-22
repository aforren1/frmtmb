# Roxygenise the packages this lane touches, into the lane library.
#   Rscript dev/shapes-roxy.R [core|sample|...]
LIB <- "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
args <- commandArgs(trailingOnly = TRUE)
which <- if (length(args)) args else c("core", "sample")
paths <- c(core = wt,
           sample = file.path(wt, "extensions/frmtmb.sample"))
for (w in which) {
  p <- if (w %in% names(paths)) paths[[w]] else file.path(wt, "extensions", w)
  cat("== roxygenise", p, "\n")
  roxygen2::roxygenise(p)   # the DESCRIPTION roclets, srr included
}
