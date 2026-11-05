# The ed() calls of review-predread-r2-release.R, on the merged source, then the install.
x <- readLines("dev/review-predread-r2-release.R")
i <- grep("^ed <- function", x)
j <- grep("^dir.create\\(lib, recursive", x)
d <- "C:/Users/adf44/source/r/rev-predread-src/merge"
lib <- "C:/Users/adf44/source/r/rev-predread-lib/merge-final"
eval(parse(text = x[i:(j - 1L)]))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_LIBS = paste(lib, "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6", sep = ";"))
inst <- function(path) {
  r <- system2(file.path(R.home("bin"), "R.exe"), c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load", paste0("--library=", lib), shQuote(path)), stdout = TRUE, stderr = TRUE)
  cat(tail(r, 1), "\n")
}
inst(d)
inst(file.path(d, "extensions/frmtmb.sample"))
