# Roxygenise and install this worktree's core (and optionally
# frmtmb.sample) into the lane's private library only.
args <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-formula2-lib"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formula2"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkgs <- if (length(args)) args else "core"
for (p in pkgs) {
  path <- if (p == "core") wt else file.path(wt, "extensions", p)
  roxygen2::roxygenise(path)
  r <- system2(file.path(R.home("bin"), "R"),
               c("CMD", "INSTALL", paste0("--library=", lib),
                 "--no-multiarch", shQuote(path)))
  if (r != 0) stop("install failed: ", p)
}
