# Roxygenise and install the lane build into the private library.
# Usage: Rscript dev/ordinal-install.R [core] [ext1 ext2 ...]
lib <- "C:/Users/adf44/source/r/wt-ordinal-lib"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
args <- commandArgs(TRUE)
if (!length(args)) args <- "core"
rcmd <- file.path(R.home("bin"), "R.exe")
for (a in args) {
  path <- if (identical(a, "core")) wt else file.path(wt, "extensions", a)
  cat("== roxygenise", path, "\n")
  roxygen2::roxygenise(path)
  cat("== install", path, "\n")
  st <- system2(rcmd, c("CMD", "INSTALL", paste0("--library=", lib),
                        "--no-multiarch", shQuote(path)))
  if (st != 0) stop("install failed: ", a)
}
cat("installed:", format(packageVersion("frmtmb", lib.loc = lib)), "\n")
