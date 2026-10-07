# Summarize dev/optima-mo-study.R output: how many fits reach the exact
# maximum, and for the ones that do not, whether the fit sits in the
# exact maximum's sign pattern (a plateau or a stop short inside the
# right cone) or in another one (another local maximum).
#   Rscript dev/optima-mo-sum.R prefix [prefix2]
args <- commandArgs(trailingOnly = TRUE)
rd <- function(p) {
  fs <- Sys.glob(paste0(p, "-*.tsv"))
  X <- do.call(rbind, lapply(fs, utils::read.delim))
  X[order(X$seed), ]
}
summ <- function(X, lab) {
  cat("==", lab, ":", nrow(X), "seeds, distinct", length(unique(X$seed)),
      "range", range(X$seed), "\n")
  g <- X$gap
  cat("  gap = exact max - logLik: <= 1e-6", sum(g <= 1e-6),
      "; > 1e-6", sum(g > 1e-6), "; > 1e-3", sum(g > 1e-3),
      "; > 1e-2", sum(g > 1e-2), "; > 0.5", sum(g > 0.5),
      "; max", format(max(g), digits = 4),
      "; min", format(min(g), digits = 4), "\n")
  cat("  optimizer code:", paste(names(table(X$code)), table(X$code),
                                 sep = "=", collapse = " "), "\n")
  same <- sign(X$b1_fit) == sign(X$b1_ex) & sign(X$b2_fit) == sign(X$b2_ex)
  miss <- g > 1e-2
  cat("  misses (> 1e-2) in the exact maximum's sign pattern:",
      sum(miss & same), "; in another pattern:", sum(miss & !same), "\n")
  cat("  min simplex weight < 1e-6 at the fit:", sum(X$min_w_fit < 1e-6),
      "; any SE not finite:", sum(X$nse_nan > 0),
      "; any warning:", sum(X$nwarn > 0), "\n")
  cat("  objective evaluations: total", sum(X$nfev), "median",
      stats::median(X$nfev), "max", max(X$nfev), "\n")
  invisible(X)
}
A <- summ(rd(args[1]), args[1])
if (length(args) > 1) {
  B <- summ(rd(args[2]), args[2])
  m <- merge(A, B, by = "seed", suffixes = c(".a", ".b"))
  cat("== paired on", nrow(m), "seeds\n")
  d <- m$ll_fit.b - m$ll_fit.a
  cat("  logLik b - a: improved > 1e-6", sum(d > 1e-6), "; worse < -1e-6",
      sum(d < -1e-6), "; max gain", format(max(d), digits = 4),
      "; max loss", format(min(d), digits = 4), "\n")
  ok_a <- m$gap.a <= 1e-6
  cat("  seeds at the maximum in a (gap <= 1e-6):", sum(ok_a),
      "; of those, |logLik b - a| max",
      format(max(abs(d[ok_a])), digits = 4), "\n")
  cat("  |exact max a - b| (identity check):",
      format(max(abs(m$ll_exact.a - m$ll_exact.b)), digits = 3), "\n")
}
