## Roxygenise and install core (and optionally extensions) into the lane's
## private library.  Rscript dev/mocode-install.R [pkg ...]
lib <- "C:/Users/adf44/source/r/wt-mocode-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-mocode"
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ";"))
pk <- commandArgs(trailingOnly = TRUE)
if (!length(pk)) pk <- "frmtmb"
r <- file.path(R.home("bin"), "R")
for (p in pk) {
  path <- if (p == "frmtmb") wt else file.path(wt, "extensions", p)
  roxygen2::roxygenise(path)
  st <- system2(r, c("CMD", "INSTALL", paste0("--library=", lib),
                     "--no-multiarch", shQuote(path)))
  if (st != 0) stop("install failed: ", p)
}
