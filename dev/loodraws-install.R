# Roxygenise and install core and/or frmtmb.sample into the lane library.
#   Rscript dev/loodraws-install.R [core] [sample]
# Give core and sample in two separate calls when core's exports change:
# roxygen loads core in this process, and the sample run then sees the
# exports of the core it loaded first.
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
lib <- "C:/Users/adf44/source/r/wt-loodraws-lib"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-loodraws"
a <- commandArgs(trailingOnly = TRUE)
R <- file.path(R.home("bin"), "R.exe")
inst <- function(path) {
  st <- system2(R, c("CMD", "INSTALL", paste0("--library=", lib),
                     "--no-multiarch", shQuote(path)))
  if (st != 0) stop("install failed: ", path)
}
if ("core" %in% a) {
  roxygen2::roxygenise(wt)
  inst(wt)
}
if ("sample" %in% a) {
  p <- file.path(wt, "extensions/frmtmb.sample")
  callr_ok <- TRUE
  roxygen2::roxygenise(p)
  inst(p)
}
