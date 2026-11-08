# Each remeasured prior pair against the record (dev/equiv-data/results.tsv,
# measured on 0.79.0): the measured verdict, the level, and the L1
# constant at the estimate.
rec <- read.delim("dev/equiv-data/results.tsv")
ids <- read.delim("dev/loglik-r80-equiv-ids.tsv")$id
rows <- lapply(ids, function(id) {
  f <- file.path("dev/loglik-equiv-out", paste0(id, ".tsv"))
  if (!file.exists(f)) {
    return(data.frame(id = id, rec_verdict = NA, now_verdict = "MISSING",
                      rec_level = NA, now_level = NA, rec_l1 = NA,
                      now_l1 = NA, same = FALSE))
  }
  now <- read.delim(f)
  r <- rec[rec$id == id, ]
  data.frame(id = id, rec_verdict = r$verdict_measured,
             now_verdict = now$verdict, rec_level = r$level,
             now_level = now$level, rec_l1 = signif(r$l1_const, 6),
             now_l1 = signif(now$l1_const, 6),
             same = identical(as.character(r$verdict_measured), as.character(now$verdict)) &&
               identical(as.character(r$level), as.character(now$level)))
})
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
cat("pairs:", nrow(tab), " verdict and level unchanged:", sum(tab$same),
    " changed:", sum(!tab$same), "\n")
