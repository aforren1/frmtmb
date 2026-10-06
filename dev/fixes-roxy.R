# Lane fixes: roxygenise core and the extensions given, against the
# lane library.
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-fixes-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in a) {
  cat("roxygenise", p, "\n")
  roxygen2::roxygenise(p)
}
