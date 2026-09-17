# Lane wt-conditions: roxygenise and install packages of this worktree
# into the lane's private library, and nowhere else.
#   Rscript dev/conditions-install.R core coupling eam ... [noroxy]
av <- commandArgs(trailingOnly = TRUE)
LIB <- "C:/Users/adf44/source/r/conditions-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
roxy <- !"noroxy" %in% av
for (p in setdiff(av, "noroxy")) {
  path <- if (p == "core") "." else file.path("extensions",
                                              paste0("frmtmb.", p))
  stopifnot(dir.exists(path))
  if (roxy) roxygen2::roxygenise(path)
  rc <- system2(file.path(R.home("bin"), "R"),
                c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
                  paste0("--library=", LIB), shQuote(path)))
  cat("R CMD INSTALL", p, "exit", rc, "\n")
  if (rc != 0) quit(status = 1)
}
