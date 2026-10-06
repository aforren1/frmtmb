# Reviewer of lane fixes, re-check: totals over dev/fixes-rev2-smooth.R
# and the fits where a build is below the best logLik of lane, base and
# mgcv ML by more than 1e-3 (a wrong answer, not a boundary).
#   Rscript dev/fixes-rev2-smooth-cmp2.R
rd <- function(lib) {
  f <- list.files("dev/fixes-rev2-log", paste0("^smooth-", lib, "-.*rds$"),
                  full.names = TRUE)
  x <- unlist(lapply(f, readRDS), recursive = FALSE)
  names(x) <- vapply(x, function(r) paste(r$case, r$seed), "")
  x
}
L <- rd("wt-fixes-lib"); B <- rd("rellib-r5")
ok <- names(L)[vapply(L, function(r) is.na(r$err), NA)]
cat("fits compared:", length(ok), "(F18 ti() refused on both)\n")
for (nm in c("lane", "base")) {
  X <- if (nm == "lane") L else B
  cat(sprintf("%s: conv != 0 %d; non-finite fixef SE %d; with a warning %d\n",
              nm, sum(vapply(X[ok], function(r) r$conv != 0, NA)),
              sum(vapply(X[ok], function(r) !r$se_ok, NA)),
              sum(vapply(X[ok], function(r) r$warns > 0, NA))))
}
worse <- list(lane = character(), base = character())
for (k in ok) {
  best <- max(L[[k]]$ll, B[[k]]$ll, L[[k]]$ml, na.rm = TRUE)
  for (nm in c("lane", "base")) {
    r <- if (nm == "lane") L[[k]] else B[[k]]
    if (best - r$ll > 1e-3) {
      worse[[nm]] <- c(worse[[nm]], sprintf("%s (%.4f, conv %d, se_ok %s)",
                                            k, best - r$ll, r$conv, r$se_ok))
    }
  }
}
for (nm in names(worse)) {
  cat(nm, "below the best by > 1e-3:", length(worse[[nm]]), "\n  ",
      paste(worse[[nm]], collapse = "\n   "), "\n")
}
