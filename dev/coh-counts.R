## The exact covered-of-total counts behind the coverage column, the
## seed range each arm used, and the selection check on the replicates
## the concurrency cap truncated. A proportion without its numerator is
## not a count, and a dropped replicate is not neutral until it is shown
## to be.
src <- readLines("dev/coh-summarize.R")
stop_at <- grep("^truth <- ", src)[1L]
eval(parse(text = paste(src[seq_len(stop_at - 1L)], collapse = "\n")))
ord <- c("cond", "smooth", "id", "idcond", "full")
for (p in c("dev/coh-recovery-main.tsv", "dev/coh-recovery-null.tsv")) {
  all <- read_tsv_kv(p)
  k <- names(which(table(all$seed) == length(unique(all$rung))))
  d <- all[as.character(all$seed) %in% k, ]
  cat("--", p, "  replicates", length(k), "  seeds", min(d$seed), "to",
      max(d$seed), "\n")
  cat("   rows", nrow(all), " every row ok:", all(all$ok),
      " partial seeds:", length(unique(all$seed)) - length(k), "\n")
  for (rg in ord) {
    s <- d[d$rung == rg, ]
    a <- all[all$rung == rg, ]
    cat(sprintf(paste0("%-8s %3d/%3d = %.4f  width %.5f",
                       "   all rows %3d/%3d = %.4f\n"),
                rg, sum(s$covers), nrow(s), mean(s$covers),
                mean(s$width), sum(a$covers), nrow(a),
                mean(a$covers)))
  }
  ## Truncation by the cap leaves a PREFIX of the fit order; attrition
  ## on hard data sets would not.
  part <- setdiff(unique(all$seed), as.numeric(k))
  for (sd_ in part) {
    got <- all$rung[all$seed == sd_]
    cat(sprintf("   partial seed %d: %s  prefix of the fit order: %s\n",
                sd_, paste(got, collapse = ","),
                identical(got, ord[seq_along(got)])))
  }
}
