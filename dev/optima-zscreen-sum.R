# Summarize dev/optima-mo-zscreen.R: per mo() term on a boundary
# simplex, the coefficient's z (marginal and conditional) on the terms
# where the full search gained over the chart-only fit, and on the rest.
X <- do.call(rbind, lapply(Sys.glob("dev/optima-log/zscreen-*.tsv"),
                           function(f) cbind(form = sub(".*zscreen-([a-z]+)-.*", "\1", f),
                                             utils::read.delim(f))))
X$gain <- X$ll_full - X$ll_chart
X$bnd <- X$min_w < 1e-6
for (fm in unique(X$form)) {
  Y <- X[X$form == fm & X$bnd, ]
  cat("==", fm, ": boundary terms", nrow(Y), "(term 1:", sum(Y$term == 1),
      ", term 2:", sum(Y$term == 2), ")\n")
  g <- Y$gain > 1e-6
  for (v in c("z_marg", "z_cond")) {
    cat(sprintf("  |%s| on gaining seeds: max %.3g | on the rest: quantiles %s\n",
                v, if (any(g)) max(abs(Y[[v]][g])) else NA,
                paste(signif(stats::quantile(abs(Y[[v]][!g]),
                                              c(0, .1, .5, .9, 1), na.rm = TRUE),
                             3), collapse = " ")))
  }
  cat("  gaining rows:", sum(g), "; seed/term/z_marg/z_cond:\n")
  if (any(g)) print(Y[g, c("seed", "term", "b", "z_marg", "z_cond", "gain")],
                    digits = 3)
}
