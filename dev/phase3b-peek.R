.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (f in commandArgs(TRUE)) {
  x <- readRDS(f)
  cat("==", f, "\n")
  print(x[setdiff(names(x), c("diagnose"))])
  cat(x$diagnose, sep = "\n")
}
