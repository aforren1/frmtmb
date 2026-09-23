d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- readRDS(paste0(d, "nomove-base.rds"))
b <- readRDS(paste0(d, "nomove-fixed.rds"))
stopifnot(identical(names(a), names(b)))
cat(sprintf("%-14s %-9s %-9s %-9s %s\n", "job", "logLik", "par", "se",
            "delta logLik"))
same <- 0L
for (nm in names(a)) {
  x <- a[[nm]]; y <- b[[nm]]
  ill <- identical(x$ll, y$ll)
  ipar <- identical(x$par, y$par)
  ise <- identical(x$se, y$se)
  if (ill && ipar && ise) same <- same + 1L
  cat(sprintf("%-14s %-9s %-9s %-9s %s\n", nm,
              if (ill) "identical" else "MOVED",
              if (ipar) "identical" else "MOVED",
              if (ise) "identical" else "MOVED",
              format(y$ll - x$ll, digits = 10)))
}
cat("\nbitwise identical on all three:", same, "of", length(a), "\n")
cat("\n-- the jobs that moved, at full precision\n")
cat(sprintf("%-14s %14s %12s %12s\n", "job", "rel d logLik",
            "max rel dpar", "max rel dse"))
rel <- function(u, v) {
  if (length(u) != length(v)) return(NA_real_)
  k <- is.finite(u) & is.finite(v)
  if (!any(k)) return(NA_real_)
  max(abs(u[k] - v[k]) / pmax(abs(u[k]), 1e-12))
}
for (nm in names(a)) {
  x <- a[[nm]]; y <- b[[nm]]
  if (identical(x$ll, y$ll) && identical(x$par, y$par) &&
        identical(x$se, y$se)) next
  cat(sprintf("%-14s %14.3e %12.3e %12.3e\n", nm,
              abs(y$ll - x$ll) / abs(x$ll), rel(x$par, y$par),
              rel(x$se, y$se)))
}
cat("jobs that moved:",
    paste(names(a)[!vapply(names(a), function(nm) {
      identical(a[[nm]]$ll, b[[nm]]$ll) && identical(a[[nm]]$par, b[[nm]]$par)
    }, TRUE)], collapse = " "), "\n")
