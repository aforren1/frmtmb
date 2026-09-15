## Roxygenise frmtmb.coupling into the worktree, loading from the
## private library first so that no other lane's build is read.
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
roxygen2::roxygenise("extensions/frmtmb.coupling")
