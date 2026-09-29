`%||%` <- function(x, y) if (is.null(x)) y else x
## Reviewer, claim 2: identical() on every stored field of the two
## identity runs. Prints a ULP-style residual where a difference exists,
## because a printed zero is not a measured zero.
a <- commandArgs(TRUE)
A <- readRDS(a[1]); B <- readRDS(a[2])
cat("rows:", length(A), length(B), "  same names:",
    identical(names(A), names(B)), "\n\n")
flds <- c("coef", "par", "logLik", "objective", "opt_obj", "convergence",
          "gmax", "grad")
allok <- TRUE
for (nm in names(A)) {
  x <- A[[nm]]; y <- B[[nm]]
  if (!is.na(x$err %||% NA) || !is.na(y$err %||% NA)) {
    cat(sprintf("%-18s ERRORED in both identically: %s\n", nm,
                identical(x$err, y$err)))
    next
  }
  bad <- character(0)
  for (f in flds) {
    if (!identical(x[[f]], y[[f]])) {
      bad <- c(bad, sprintf("%s(maxdiff=%s)", f,
                            format(max(abs(unlist(x[[f]]) -
                                             unlist(y[[f]]))), digits = 6)))
    }
  }
  wsame <- identical(x$warns, y$warns)
  cat(sprintf("%-18s fit bitwise identical: %-5s  warnings identical: %-5s%s\n",
              nm, length(bad) == 0, wsame,
              if (length(bad)) paste0("  DIFF: ", paste(bad, collapse = " "))
              else ""))
  if (length(bad)) allok <- FALSE
  if (!wsame) {
    cat("   base-only warnings:\n")
    for (w in setdiff(y$warns, x$warns)) cat("     -", w, "\n")
    cat("   lane-only warnings:\n")
    for (w in setdiff(x$warns, y$warns)) cat("     +", w, "\n")
  }
}
cat("\nALL FITS BITWISE IDENTICAL:", allok, "\n")
