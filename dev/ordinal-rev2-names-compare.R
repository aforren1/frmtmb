# Reviewer re-check: identical() of every exposed name, base vs lane.
# Output: dev/ordinal-rev2-log-names-compare.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
b <- readRDS(file.path(wt, "dev/ordinal-rev2-names-base.rds"))
l <- readRDS(file.path(wt, "dev/ordinal-rev2-names-lane.rds"))
sink(file.path(wt, "dev/ordinal-rev2-log-names-compare.txt"), split = TRUE)
tot <- 0; same <- 0
for (m in names(b)) {
  ks <- union(names(b[[m]]), names(l[[m]]))
  s <- vapply(ks, function(k) identical(b[[m]][[k]], l[[m]][[k]]), NA)
  errs <- ks[vapply(ks, function(k) {
    x <- l[[m]][[k]]; is.character(x) && length(x) == 1 && startsWith(x, "ERROR")
  }, NA)]
  tot <- tot + length(ks); same <- same + sum(s)
  cat(sprintf("%-18s %2d/%2d identical  differ: %s  errors(both arms): %s\n", m,
              sum(s), length(ks), paste(ks[!s], collapse = ","),
              paste(errs, collapse = ",")))
  zl <- l[[m]]$template_lengths
  if (is.numeric(zl) && any(zl == 0)) cat("   zero-length template blocks:",
                                          names(zl)[zl == 0], "\n")
}
cat("\ntotal:", same, "of", tot, "identical\n")
for (m in names(b)) for (k in union(names(b[[m]]), names(l[[m]])))
  if (!identical(b[[m]][[k]], l[[m]][[k]])) {
    cat("\n----", m, k, "\nbase:\n"); print(utils::head(b[[m]][[k]], 12))
    cat("lane:\n"); print(utils::head(l[[m]][[k]], 12))
  }
sink()
