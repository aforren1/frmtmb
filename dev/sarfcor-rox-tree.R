# Roxygenise core, frmtmb.spatial and frmtmb.sample in a tree outside
# the worktree (the part-A stage tree of dev/sarfcor-ab-diff.sh).
#   Rscript dev/sarfcor-rox-tree.R <tree>
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
tree <- commandArgs(TRUE)[1]
for (d in c("", "extensions/frmtmb.spatial", "extensions/frmtmb.sample")) {
  p <- normalizePath(file.path(tree, d), winslash = "/")
  cat("roxygenise", p, "\n")
  suppressMessages(roxygen2::roxygenise(p))
}
cat("ROX DONE\n")
