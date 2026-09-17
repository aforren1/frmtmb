## Compare dev/famlink-rev2-silent-{base,lane}.rds.
b <- readRDS("dev/famlink-rev2-silent-base.rds"); l <- readRDS("dev/famlink-rev2-silent-lane.rds")
for (k in names(l)) {
  f <- c("ll", "est", "se", "conv", "msg", "error", "vcov_warnings", "diag", "diag_warn", "print", "summary")
  same <- vapply(f, function(z) identical(b[[k]][[z]], l[[k]][[z]]), TRUE)
  nw_b <- sum(grepl("NA/NaN", b[[k]]$warnings)); nw_l <- sum(grepl("NA/NaN", l[[k]]$warnings))
  cat(sprintf("%-24s differs in: %-40s NaN warnings base %2d lane %2d; lane fit$opt$nonfinite_trials %s; other warnings identical %s\n",
      k, paste(f[!same], collapse = ","), nw_b, nw_l, format(l[[k]]$nonfinite),
      identical(b[[k]]$warnings[!grepl("NA/NaN", b[[k]]$warnings)], l[[k]]$warnings[!grepl("NA/NaN", l[[k]]$warnings)])))
}
