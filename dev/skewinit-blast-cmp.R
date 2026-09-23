d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- readRDS(paste0(d, "blast-base.rds"))
b <- readRDS(paste0(d, "blast-fixed.rds"))
stopifnot(identical(names(a), names(b)))
rel <- function(u, v) {
  if (length(u) != length(v)) return(NA_real_)
  k <- is.finite(u) & is.finite(v)
  if (!any(k)) return(NA_real_)
  max(abs(u[k] - v[k]) / pmax(abs(u[k]), 1e-12))
}
cat(sprintf("%-16s %-9s %12s %12s %12s %s\n", "job", "bitwise",
            "rel dlogLik", "rel dpar", "rel dse", "conv/warn"))
same <- 0L
for (nm in names(a)) {
  x <- a[[nm]]; y <- b[[nm]]
  bit <- identical(x$ll, y$ll) && identical(x$par, y$par) &&
    identical(x$se, y$se)
  if (bit) same <- same + 1L
  cat(sprintf("%-16s %-9s %12.3e %12.3e %12.3e %d/%d -> %d/%d\n", nm,
              if (bit) "identical" else "moved",
              abs(y$ll - x$ll) / max(abs(x$ll), 1e-12),
              rel(x$par, y$par), rel(x$se, y$se),
              x$conv, x$warns, y$conv, y$warns))
}
cat("\nbitwise identical:", same, "of", length(a), "\n")
worse <- vapply(names(a), function(nm) b[[nm]]$ll < a[[nm]]$ll - 1e-8, TRUE)
cat("jobs where the lane is LOWER by more than 1e-8:", sum(worse),
    if (any(worse)) paste(names(a)[worse], collapse = " "), "\n")
better <- vapply(names(a), function(nm) b[[nm]]$ll > a[[nm]]$ll + 1e-8, TRUE)
cat("jobs where the lane is HIGHER by more than 1e-8:", sum(better),
    if (any(better)) paste(names(a)[better], collapse = " "), "\n")
for (nm in names(a)[better]) {
  cat(sprintf("   %s: base %.9f  lane %.9f  gain %.6f\n", nm,
              a[[nm]]$ll, b[[nm]]$ll, b[[nm]]$ll - a[[nm]]$ll))
}
cat("new warnings introduced:",
    sum(vapply(names(a), function(nm) b[[nm]]$warns > a[[nm]]$warns, TRUE)),
    "\n")
cat("new nonzero convergence codes:",
    sum(vapply(names(a),
               function(nm) b[[nm]]$conv != 0 & a[[nm]]$conv == 0, TRUE)),
    "\n")
