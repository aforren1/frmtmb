# Roxygenise one package of the worktree with the lane's library first.
#   Rscript dev/mvnlatent-roxy.R <package dir>
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
roxygen2::roxygenise(commandArgs(trailingOnly = TRUE)[1])
