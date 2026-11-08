# Fix round M1: pkgdown::check_pkgdown() on the worktree (reads only).
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("lib: pkgdown", format(packageVersion("pkgdown")),
    dirname(find.package("pkgdown")), "\n")
r <- tryCatch({
  pkgdown::check_pkgdown("C:/Users/adf44/source/r/frmtmb-wt-loglik")
  "OK"
}, error = function(e) conditionMessage(e))
cat("check_pkgdown:", r, "\n")
