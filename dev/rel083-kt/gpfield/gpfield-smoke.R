# Lane gpfield: each row's default fit, with its warnings' classes, theta
# and the second descents' count. Args: arm tags(comma)
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/gpfield/gpfield-lib.R")
gf_libs(a[1])
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/gpfield/gpfield-rows.R")
for (tag in strsplit(a[2], ",")[[1]]) {
  r <- gf_row(tag)
  t0 <- proc.time()[[3]]
  q <- gf_fit(r$fml, data = r$data, family = r$family)
  wall <- proc.time()[[3]] - t0
  f <- q$fit
  if (inherits(f, "error")) {
    cat(sprintf("SMOKE %s ERROR %s\n", tag, conditionMessage(f)))
    next
  }
  gs <- f$opt$gp_second
  cat(sprintf("SMOKE %s ll %.6f conv %d nwarn %d wall %.2f runs %s gain %s theta %s | %s\n",
              tag, as.numeric(logLik(f)), f$opt$convergence, length(q$warn),
              wall, if (is.null(gs)) "NA" else gs[["runs"]],
              if (is.null(gs)) "NA" else format(signif(gs[["gain"]], 4)),
              paste(sprintf("%.3f", f$estimates$theta), collapse = ","),
              paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
  for (w in q$warn) cat("  W", substr(w, 1, 300), "\n")
  flush.console()
}
cat("DONE\n")
