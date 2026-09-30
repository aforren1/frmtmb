# Reviewer, lane ordinal: identical() of every output of
# dev/ordinal-rev-bitwise.R, base against lane.
# Output: dev/ordinal-rev-log-bitwise-compare.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
b <- readRDS(file.path(wt, "dev/ordinal-rev-bitwise-base.rds"))
l <- readRDS(file.path(wt, "dev/ordinal-rev-bitwise-lane.rds"))
sink(file.path(wt, "dev/ordinal-rev-log-bitwise-compare.txt"), split = TRUE)
cat("models:", length(b), "base,", length(l), "lane\n")
tab <- NULL
for (m in names(b)) {
  keys <- union(names(b[[m]]), names(l[[m]]))
  same <- vapply(keys, function(k) identical(b[[m]][[k]], l[[m]][[k]]), NA)
  err <- vapply(keys, function(k) {
    x <- l[[m]][[k]]
    is.character(x) && length(x) == 1L && startsWith(x, "ERROR")
  }, NA)
  tab <- rbind(tab, data.frame(model = m, n = length(keys),
                               identical = sum(same),
                               differ = paste(keys[!same], collapse = ","),
                               lane_errors = paste(keys[err], collapse = ",")))
}
print(tab, row.names = FALSE, right = FALSE)
cat("\noutputs compared:", sum(tab$n), " identical:", sum(tab$identical), "\n")
# what differs, shown
show <- function(m, k) {
  cat("\n----", m, k, "\n")
  x <- b[[m]][[k]]; y <- l[[m]][[k]]
  if (is.character(x) && is.character(y)) {
    cat("base only:\n"); print(setdiff(x, y))
    cat("lane only:\n"); print(setdiff(y, x))
  } else {
    cat("base:\n"); print(utils::head(x, 8))
    cat("lane:\n"); print(utils::head(y, 8))
  }
}
for (i in seq_len(nrow(tab))) {
  ks <- strsplit(tab$differ[i], ",")[[1]]
  for (k in ks) show(tab$model[i], k)
}
sink()
