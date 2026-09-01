## Install the wt-benchp worktree copy of frmtmb into a scratch library,
## used only to obtain the reference frm() InstEval logLik.
scratch <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/libp"
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
install.packages("C:/Users/adf44/source/r/frmtmb-wt-benchp", repos = NULL,
                 type = "source", lib = scratch, INSTALL_opts = "--no-docs")
.libPaths(c(scratch, .libPaths()))
cat("frmtmb:", as.character(utils::packageVersion("frmtmb")), "\n")
