## Roxygenise, then install into the LANE-PRIVATE library only.
## Usage: Rscript dev/argspell-install.R [pkgdir ...]
LIB <- "C:/Users/adf44/source/r/argspell-lib"
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- "."
for (d in args) {
  cat("\n==== roxygenise", d, "====\n")
  roxygen2::roxygenise(d)
  cat("==== install", d, "====\n")
  r <- system2(file.path(R.home("bin"), "R"),
               c("CMD", "INSTALL", paste0("--library=", shQuote(LIB)),
                 "--no-multiarch", shQuote(normalizePath(d))),
               stdout = TRUE, stderr = TRUE)
  cat(tail(r, 12), sep = "\n")
  if (!is.null(attr(r, "status")) && attr(r, "status") != 0) {
    cat(r, sep = "\n")
    stop("install failed: ", d)
  }
}
cat("\nlibrary now holds:\n")
print(utils::installed.packages(lib.loc = LIB)[, "Version", drop = FALSE])
