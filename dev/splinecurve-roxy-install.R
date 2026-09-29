# Lane splinecurve: roxygenise frmtmb.spline and install it into the
# private library, never anywhere else.
lib <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- "C:/Users/adf44/source/r/frmtmb-wt-release/extensions/frmtmb.spline"
stopifnot(packageVersion("frmtmb", lib.loc = lib) >= "0.65.0")
roxygen2::roxygenise(pkg)
st <- system2(file.path(R.home("bin"), "R.exe"),
              c("CMD", "INSTALL", paste0("--library=", lib), "--no-multiarch",
                shQuote(pkg)))
if (st != 0) stop("install failed with status ", st)
cat("installed frmtmb.spline",
    format(packageVersion("frmtmb.spline", lib.loc = lib)), "into", lib, "\n")
