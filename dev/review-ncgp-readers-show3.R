source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
a <- readRDS(file.path(REV_OUT, "readers-lane.rds")); b <- readRDS(file.path(REV_OUT, "readers-base.rds"))
rel <- function(x, y) { x <- as.numeric(unlist(x)); y <- as.numeric(unlist(y)); ok <- is.finite(x) & is.finite(y); if (length(x) != length(y) || !any(ok)) return(NA); m <- which.max(abs(x[ok] - y[ok])); c(abs = abs(x[ok] - y[ok])[m], at = y[ok][m]) }
for (k in c("deriv", "deriv2", "emtr", "fitted", "fittedn", "linp", "ran", "resid", "ci", "se", "curve", "mpred")) {
  cat("==", k, "\n")
  for (m in names(a)) {
    x <- a[[m]][[k]]; y <- b[[m]][[k]]
    if (is.character(x) || is.character(y)) next
    if (k %in% c("deriv", "deriv2")) { x <- x[, c(".estimate", ".se")]; y <- y[, c(".estimate", ".se")] }
    r <- rel(x, y)
    cat(sprintf("  %-8s max abs diff %.3e where base reads %.4g\n", m, r[1], r[2]))
  }
}
