# Reviewer: re-roxygenise a copy of the worktree, to check that man/ and
# NAMESPACE in the lane's tree are current.
.libPaths(c("C:/Users/adf44/source/r/cifixrev-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- "C:/Users/adf44/source/r/frmtmb-wt-cifix/dev/cifixrev-check/src"
roxygen2::roxygenise(src)
roxygen2::roxygenise(file.path(src, "extensions/frmtmb.spline"))
