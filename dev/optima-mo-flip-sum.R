# Summarize dev/optima-mo-flip.R: seeds reaching the exact maximum
# (gap <= 1e-6) per arm, and the extra objective evaluations.
#   Rscript dev/optima-mo-flip-sum.R prefix
args <- commandArgs(trailingOnly = TRUE)
X <- do.call(rbind, lapply(Sys.glob(paste0(args[1], "-*.tsv")),
                           utils::read.delim))
cat("==", args[1], ":", nrow(X), "seeds, distinct", length(unique(X$seed)),
    "\n")
for (v in c("gap_fit", "gap_flip", "gap_vertex", "gap_both")) {
  g <- X[[v]]
  cat(sprintf("  %-10s <= 1e-6 %3d | > 1e-6 %3d | > 1e-2 %3d | > 0.5 %3d | max %.4g\n",
              v, sum(g <= 1e-6), sum(g > 1e-6), sum(g > 1e-2),
              sum(g > 0.5), max(g)))
}
cat("  evaluations: fit total", sum(X$ev_fit), "| flip extra",
    sum(X$ev_flip), "| vertex extra", sum(X$ev_vertex), "| both extra",
    sum(X$ev_both), "\n")
cat("  autoscaled fits:", sum(X$autoscaled), "\n")
