# Lane ndaterm: roxygenise and install core and the changed extensions
# into the lane library. Usage: Rscript dev/ndaterm-install.R [pkgs...]
lib <- "C:/Users/adf44/source/r/wt-ndaterm-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ndaterm"
pk <- commandArgs(TRUE)
if (!length(pk)) pk <- "frmtmb"
for (p in pk) {
  dir <- if (p == "frmtmb") wt else file.path(wt, "extensions", p)
  roxygen2::roxygenise(dir)
  rc <- system2(file.path(R.home("bin"), "R"),
                c("CMD", "INSTALL", paste0("--library=", lib),
                  "--no-multiarch", shQuote(dir)))
  if (rc != 0) stop("install of ", p, " failed")
  cat("installed", p, "\n")
}
