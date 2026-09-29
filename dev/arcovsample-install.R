# Lane wt-arcovsample: roxygenise then install core and frmtmb.sample
# into the LANE-PRIVATE library. Nothing else is written to.
#
#   Rscript dev/arcovsample-install.R

LIB <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
WT <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample"
dir.create(LIB, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(LIB, USER))

stopifnot(identical(normalizePath(.libPaths()[1], "/"),
                    normalizePath(LIB, "/")))

roxygen2::roxygenise(WT)
roxygen2::roxygenise(file.path(WT, "extensions/frmtmb.sample"))

R <- file.path(R.home("bin"), "R")
for (p in c(WT, file.path(WT, "extensions/frmtmb.sample"))) {
  st <- system2(R, c("CMD", "INSTALL", paste0("--library=", LIB),
                     "--no-multiarch", shQuote(p)))
  cat("INSTALL", basename(p), "status", st, "\n")
  if (st != 0L) stop("install failed: ", p)
}
cat("frmtmb        ", format(packageVersion("frmtmb", lib.loc = LIB)), "\n")
cat("frmtmb.sample ",
    format(packageVersion("frmtmb.sample", lib.loc = LIB)), "\n")
