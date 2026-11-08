a <- commandArgs(TRUE)
p <- function(t) readRDS(sprintf("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-noreg-%s.rds", t))
A <- p(a[1]); B <- p(a[2]); C <- if (length(a) > 2) p(a[3])
cat(sprintf("%-12s %-9s %-13s %-12s %-10s %s\n", "design", "basis", "par identical",
            "max|dpar|rel", "dobj", "SE max rel"))
nid <- 0
for (nm in names(A)) {
  x <- A[[nm]]; y <- B[[nm]]
  if (!is.null(x$err) || !is.null(y$err)) { cat(nm, "error\n"); next }
  idp <- identical(unname(x$par), unname(y$par)) && identical(x$obj, y$obj)
  nid <- nid + idp
  dp <- max(abs(x$par - y$par) / pmax(abs(y$par), 1e-8))
  ds <- max(abs(x$se - y$se) / pmax(abs(y$se), 1e-300))
  ctrl <- if (!is.null(C)) identical(unname(B[[nm]]$par), unname(C[[nm]]$par)) else NA
  cat(sprintf("%-12s %-9s %-13s %-12.3g %-10.3g %-10.3g base/base identical %s\n", nm,
              x$basis, idp, dp, x$obj - y$obj, ds, ctrl))
}
cat("identical", nid, "of", length(A), "\n")
