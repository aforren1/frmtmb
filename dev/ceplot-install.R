# Lane ceplot: roxygenise and install core, then frmtmb.sample, into the
# lane's private library.
#   Rscript dev/ceplot-install.R [core|sample|both]
what <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(what)) what <- "both"
LIB <- "C:/Users/adf44/source/r/wt-ceplot-lib"
WT <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot"
dir.create(LIB, showWarnings = FALSE)
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(PATH = paste("C:/rtools45/usr/bin",
                        "C:/rtools45/x86_64-w64-mingw32.static.posix/bin",
                        Sys.getenv("PATH"), sep = ";"))
rcmd <- file.path(R.home("bin"), "R")
inst <- function(path) {
  st <- system2(rcmd, c("CMD", "INSTALL", paste0("--library=", LIB),
                        "--no-multiarch", shQuote(path)))
  if (st != 0) stop("install failed: ", path)
}
if (what %in% c("core", "both")) {
  roxygen2::roxygenise(WT)
  inst(WT)
}
if (what %in% c("sample", "both")) {
  # the sample package's roxygen loads frmtmb from the library, so core
  # goes in first
  roxygen2::roxygenise(file.path(WT, "extensions", "frmtmb.sample"))
  inst(file.path(WT, "extensions", "frmtmb.sample"))
}
cat("installed", what, "into", LIB, "\n")
