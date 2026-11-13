# Lane ci082: roxygenise the worktree's core, then install it into the
# lane library only.
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ci082"
roxygen2::roxygenise(wt)
r <- system2(file.path(R.home("bin"), "R.exe"),
             c("CMD", "INSTALL", "--library=C:/Users/adf44/source/r/wt-ci082-lib",
               "--no-multiarch", shQuote(wt)))
cat("INSTALL status", r, "\n")
