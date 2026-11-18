# Lane gpselimit: bitwise comparison of the sweep's saved rows between two
# arms (dev/gpselimit-out/sweep-<a>-<row>.rds against sweep-<b>-<row>.rds):
# estimates, outer point, logLik, every SE (NaN matched as NaN), the lost
# reasons, and the conditions (classes).
# Usage: Rscript gpselimit-cmp.R <a> <b> [tag]
a <- commandArgs(TRUE)
tag <- if (length(a) >= 3) a[3] else ""
out <- "C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-out"
same <- function(x, y) identical(x, y)
for (row in 1:30) {
  fa <- file.path(out, sprintf("sweep-%s%s-%02d.rds", a[1], tag, row))
  fb <- file.path(out, sprintf("sweep-%s%s-%02d.rds", a[2], tag, row))
  if (!file.exists(fa) || !file.exists(fb)) {
    cat("CMP", row, "missing\n"); next
  }
  x <- readRDS(fa); y <- readRDS(fb)
  cat(sprintf(paste0("CMP %2d est %d par %d ll %d se %d lost %d ",
                     "wclass %d mclass %d | %s -> %s | %s -> %s\n"),
              row, same(x$est, y$est), same(x$par, y$par), same(x$ll, y$ll),
              same(x$se, y$se), same(x$lost, y$lost),
              same(x$wclass, y$wclass), same(x$mclass, y$mclass),
              paste(unique(x$wclass), collapse = ";"),
              paste(unique(y$wclass), collapse = ";"),
              paste(unique(x$mclass), collapse = ";"),
              paste(unique(y$mclass), collapse = ";")))
}
