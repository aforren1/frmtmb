# dev/machine-library.md: verify a restored library with a FIT, not a
# version string. The reference is this lane's own 27-job base battery,
# measured at 15:14 against rellib-r3 BEFORE the 17:05 loss.
d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- readRDS(paste0(d, "nomove-base.rds"))
b <- readRDS(paste0(d, "nomove-base-postrestore.rds"))
stopifnot(identical(names(a), names(b)))
same <- 0L
moved <- character(0)
for (nm in names(a)) {
  if (identical(a[[nm]]$ll, b[[nm]]$ll) &&
        identical(a[[nm]]$par, b[[nm]]$par) &&
        identical(a[[nm]]$se, b[[nm]]$se)) {
    same <- same + 1L
  } else {
    moved <- c(moved, nm)
  }
}
cat("jobs bitwise identical before and after the restore:", same, "of",
    length(a), "\n")
if (length(moved)) {
  cat("MOVED:", paste(moved, collapse = " "), "\n")
  for (nm in moved) {
    cat(sprintf("  %-14s d logLik %s\n", nm,
                format(b[[nm]]$ll - a[[nm]]$ll, digits = 10)))
  }
} else {
  cat("the restored toolchain reproduces every pre-loss number exactly\n")
}
op <- options(digits = 12)
cat("\ntwo reference values, for the next restore:\n")
for (nm in c("gaussian_re", "skew_icpt")) {
  cat(sprintf("  %-14s logLik %s\n", nm, format(b[[nm]]$ll)))
}
options(op)
