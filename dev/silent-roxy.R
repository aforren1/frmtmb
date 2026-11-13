# Roxygenises one package of the worktree with the lane library first.
# Usage: Rscript dev/silent-roxy.R <package dir>
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
roxygen2::roxygenise(commandArgs(TRUE)[1])
