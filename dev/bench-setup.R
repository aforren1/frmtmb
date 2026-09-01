## one-time setup: install optimParallel (CRAN) and the worktree copy of
## frmtmb into a scratch library so PSOCK workers can library() it.
scratch <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib"
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("optimParallel", quietly = TRUE)) {
  install.packages("optimParallel", repos = "https://cloud.r-project.org")
}
cat("optimParallel:", as.character(utils::packageVersion("optimParallel")),
    "\n")

install.packages("C:/Users/adf44/source/r/frmtmb-wt-bench", repos = NULL,
                 type = "source", lib = scratch, INSTALL_opts = "--no-docs")
.libPaths(c(scratch, .libPaths()))
cat("frmtmb:", as.character(utils::packageVersion("frmtmb")), "\n")
cat("frmtmb lib:", dirname(system.file(package = "frmtmb")), "\n")
