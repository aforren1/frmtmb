# Roxygenise core and install it into the lane's private library ONLY.
# Run: Rscript dev/famlink-install.R [extension ...]
lib <- "C:/Users/adf44/source/r/famlink-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(identical(normalizePath(.libPaths()[1], winslash = "/"),
                    normalizePath(lib, winslash = "/")))
root <- "C:/Users/adf44/source/r/frmtmb-wt-famlink"
pkgs <- commandArgs(trailingOnly = TRUE)
if (!length(pkgs)) pkgs <- "frmtmb"
for (p in pkgs) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  roxygen2::roxygenise(path)
  st <- system2(file.path(R.home("bin"), "R.exe"),
                c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
                  paste0("--library=", shQuote(lib)), shQuote(path)))
  if (st != 0) stop("install of ", p, " failed")
}
cat("installed", pkgs, "into", lib, "\n")
