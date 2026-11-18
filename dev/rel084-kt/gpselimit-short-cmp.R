# rel084-kt copy of dev/reviews/gpselimit-short-cmp.R: b = argument 1 (an a-*.rds path), l = argument 2; the lane copy read a-base.rds and a-lane.rds from the lane tree. No other change.
# Short check of lane gpselimit's fix round 1: base against lane, field by
# field, on gpselimit-short-a.R's fits (separate processes).
# Usage: Rscript gpselimit-short-cmp.R
out <- "C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-short-out"
b <- readRDS(commandArgs(TRUE)[1])
l <- readRDS(commandArgs(TRUE)[2])
for (n in names(b)) {
  x <- b[[n]]; y <- l[[n]]
  cat(sprintf("CMP %-10s resolved %d | est %d par %d ll %d se %d lost %d\n",
              n, !any(x$white), identical(x$est, y$est),
              identical(x$par, y$par), identical(x$ll, y$ll),
              identical(unname(x$se), unname(y$se)),
              identical(unname(x$lost), unname(y$lost))))
}
