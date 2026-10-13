# Reviewer of lane seflags, final check: compare dev/seflags-rev3-bitwise.R
# outputs of two arms, element by element, with identical().
#   Rscript dev/seflags-rev3-bitcmp.R <a.rds> <b.rds>
a <- readRDS(commandArgs(TRUE)[1]); b <- readRDS(commandArgs(TRUE)[2])
for (k in names(a)) {
  x <- a[[k]]; y <- b[[k]]
  same <- vapply(c("fixef", "vcov", "fitted", "fitted_pop", "predict_lin"),
                 function(e) identical(x[[e]], y[[e]]), NA)
  rel <- vapply(c("fixef", "fitted", "fitted_pop", "predict_lin"), function(e) {
    if (is.null(x[[e]]) || is.null(y[[e]])) return(NA_real_)
    max(abs(x[[e]] / y[[e]] - 1), na.rm = TRUE)
  }, 0)
  cat(sprintf("%-10s lost %d/%d repaired %s/%s | identical %s | max rel diff %s\n",
              k, length(x$lost), length(y$lost), x$repaired, y$repaired,
              if (all(same)) "all" else paste(names(same)[!same], collapse = ","),
              paste(signif(rel, 2), collapse = " ")))
}
