# Roxygenise and install this lane's packages into the private library.
# Usage: Rscript dev/postfit2-install.R [core|sample|both]
what <- commandArgs(TRUE)[1]
if (is.na(what)) what <- "both"
lib <- "C:/Users/adf44/source/r/wt-postfit2-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
rcmd <- file.path(R.home("bin"), "R.exe")
inst <- function(path) {
  st <- system2(rcmd, c("CMD", "INSTALL", paste0("--library=", lib),
                        "--no-multiarch", shQuote(path)))
  if (st != 0) stop("install failed: ", path)
}
if (what %in% c("core", "both")) {
  roxygen2::roxygenise(wt)
  inst(wt)
}
if (what %in% c("sample", "both")) {
  roxygen2::roxygenise(file.path(wt, "extensions/frmtmb.sample"))
  inst(file.path(wt, "extensions/frmtmb.sample"))
}
