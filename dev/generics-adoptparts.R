# Where do the 8.89 us of the adopted path actually go?  Reported
# because the number is large enough to want an attribution rather
# than a shrug, and because the cheap-looking call may not be cheap.
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(loadNamespace("loo"))
lns <- asNamespace("loo")
memo <- lns
fb <- function(x, ...) UseMethod("loo")

block <- function(fun, reps) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) fun()
  proc.time()[["elapsed"]] - t0
}
grow <- function(fun, label) {
  reps <- 1024L
  repeat {
    el <- block(fun, reps)
    if (el >= 1.2 || reps > 8e6) break
    reps <- reps * 4L
  }
  best <- Inf
  for (r in 1:5) best <- min(best, block(fun, reps))
  cat(sprintf("%-40s %9d reps %7.3f s %8.2f us\n", label, reps, best,
              best / reps * 1e6))
  invisible(best / reps)
}
cat("```\n")
cat("== where the adopted path's time goes ==\n")
grow(function() isNamespaceLoaded("loo"), "isNamespaceLoaded('loo')")
grow(function() asNamespace("loo"), "asNamespace('loo')")
grow(function() identical(lns, memo), "identical(ns, memo)")
grow(function() getExportedValue("loo", "loo"), "getExportedValue()")
grow(function() {
  if (isNamespaceLoaded("loo")) {
    ons <- asNamespace("loo")
    if (!identical(ons, memo)) memo <<- ons
    fb
  } else fb
}, "the whole accessor body, inline")
cat("```\n")
