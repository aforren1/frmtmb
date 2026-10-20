# Lane adhess: compare dev/adhess-move.R's two arms fit by fit.
#   Rscript dev/adhess-movecmp.R <lane.rds> <base.rds>
a <- commandArgs(TRUE)
L <- readRDS(a[1])
B <- readRDS(a[2])
cov_err <- function(V, Vr) {
  dg <- sqrt(diag(Vr))
  c(se = max(abs(sqrt(diag(V)) / dg - 1)),
    cor = max(abs(V - Vr) / outer(dg, dg)))
}
rows <- list()
for (k in names(L)) {
  l <- L[[k]]; b <- B[[k]]
  if (!is.null(l$error) || !is.null(b$error) || is.null(l$V) ||
        is.null(b$V)) {
    cat(sprintf("%-12s skipped (error or no vcov)\n", k))
    next
  }
  same_par <- identical(unname(l$par), unname(b$par))
  mv <- if (all(is.finite(l$V)) && all(is.finite(b$V))) cov_err(l$V, b$V) else
    c(se = NA, cor = NA)
  el <- if (!is.null(l$ref)) cov_err(l$V, l$ref) else c(se = NA, cor = NA)
  eb <- if (!is.null(l$ref)) cov_err(b$V, l$ref) else c(se = NA, cor = NA)
  rows[[k]] <- data.frame(fit = k, re = l$re, same_par = same_par,
                          move_se = mv[["se"]], move_cor = mv[["cor"]],
                          base_se = eb[["se"]], lane_se = el[["se"]],
                          base_cor = eb[["cor"]], lane_cor = el[["cor"]],
                          warn_same = identical(l$warn, b$warn),
                          lost_same = identical(l$lost, b$lost),
                          V_same = identical(l$V, b$V))
}
tab <- do.call(rbind, rows)
op <- options(width = 200)
print(format(tab, digits = 3), row.names = FALSE)
for (re in c(FALSE, TRUE)) {
  t1 <- tab[tab$re == re, ]
  cat(sprintf(paste0("\n%s random effects: %d fits | estimates identical %d | ",
                     "warnings identical %d | lost sets identical %d\n",
                     "  move SE: max %.2e median %.2e | move cor: max %.2e\n",
                     "  vs reference SE: base max %.2e median %.2e | ",
                     "lane max %.2e median %.2e\n",
                     "  lane closer on %d of %d | vcov bitwise identical %d\n"),
              if (re) "WITH" else "WITHOUT", nrow(t1), sum(t1$same_par),
              sum(t1$warn_same), sum(t1$lost_same),
              max(t1$move_se, na.rm = TRUE), stats::median(t1$move_se, na.rm = TRUE),
              max(t1$move_cor, na.rm = TRUE),
              max(t1$base_se, na.rm = TRUE), stats::median(t1$base_se, na.rm = TRUE),
              max(t1$lane_se, na.rm = TRUE), stats::median(t1$lane_se, na.rm = TRUE),
              sum(t1$lane_se < t1$base_se, na.rm = TRUE),
              sum(is.finite(t1$lane_se) & is.finite(t1$base_se)),
              sum(t1$V_same)))
}
options(op)
