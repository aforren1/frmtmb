# Install the worktree's core (and named extensions) into the lane's
# private library: roxygenise, then R CMD INSTALL --library=LIB.
# Usage: Rscript dev/ciharden-install.R [core] [ext ...]
LIB <- "C:/Users/adf44/source/r/wt-ciharden-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-ciharden"
a <- commandArgs(TRUE)
for (p in a) {
  dir <- if (p == "core") root else file.path(root, "extensions", p)
  roxygen2::roxygenise(dir)
  r <- system2(file.path(R.home("bin"), "R.exe"),
               c("CMD", "INSTALL", paste0("--library=", LIB),
                 "--no-multiarch", shQuote(dir)))
  if (r != 0) stop("install failed: ", p)
}
