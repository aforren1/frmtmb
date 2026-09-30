# Lane sampfix: roxygenise and install core and frmtmb.sample into the
# private library.   Rscript dev/sampfix-install.R [core] [sample]
a <- commandArgs(trailingOnly = TRUE)
LIB <- "C:/Users/adf44/source/r/wt-sampfix-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix"
R <- file.path(R.home("bin"), "R.exe")
inst <- function(path) {
  st <- system2(R, c("CMD", "INSTALL", paste0("--library=", LIB),
                     "--no-multiarch", shQuote(path)),
                stdout = TRUE, stderr = TRUE)
  cat(tail(st, 2), sep = "\n")
}
if ("core" %in% a) {
  roxygen2::roxygenise(WT)
  inst(WT)
}
if ("sample" %in% a) {
  roxygen2::roxygenise(file.path(WT, "extensions/frmtmb.sample"))
  inst(file.path(WT, "extensions/frmtmb.sample"))
}
