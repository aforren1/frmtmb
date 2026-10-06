# Reviewer round 2: roxygenise and install the trial merge of lanes gpby
# and fixes (scratchpad tree) into the reviewer's own library.
tree <- commandArgs(TRUE)[1]
lib <- "C:/Users/adf44/source/r/gpby-rev-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
roxygen2::roxygenise(tree)
R <- file.path(R.home("bin"), "R.exe")
for (p in c(".", "extensions/frmtmb.spline", "extensions/frmtmb.sample")) {
  if (p != ".") roxygen2::roxygenise(file.path(tree, p))
  st <- system2(R, c("CMD", "INSTALL", paste0("--library=", lib),
                     "--no-multiarch", shQuote(file.path(tree, p))))
  cat("INSTALL", p, "status", st, "\n")
}
