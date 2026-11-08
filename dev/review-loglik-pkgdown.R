# Review of lane loglik: does the pkgdown reference index list every
# topic? Runs pkgdown::check_pkgdown() (reads only) on the worktree and
# on main at 3cd2a0c8 for comparison.
.libPaths(c("C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("lib: pkgdown", format(packageVersion("pkgdown")),
    dirname(find.package("pkgdown")), "\n")
for (p in c("C:/Users/adf44/source/r/frmtmb-wt-loglik",
            "C:/Users/adf44/source/r/frmtmb")) {
  cat("==", p, "\n")
  r <- tryCatch({ pkgdown::check_pkgdown(p); "OK" },
                error = function(e) conditionMessage(e))
  cat(r, "\n")
}
