# Every sampled pair's live frm_compat() status beside the status frozen
# in dev/equiv-data/sample.tsv (its stratum), with no fit. Prints the pairs
# where they differ, with the record's verdict for each.
#   Rscript dev/compat-live-status.R <lib;lib> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
smp <- read.delim("dev/equiv-data/sample.tsv", stringsAsFactors = FALSE)
rec <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
cp <- as.data.frame(frm_compat())
live <- vapply(seq_len(nrow(smp)), function(i) {
  a <- smp$feature_a[i]
  b <- smp$feature_b[i]
  s <- cp$status[(cp$feature_a == a & cp$feature_b == b) |
                   (cp$feature_a == b & cp$feature_b == a)]
  if (length(s) != 1L) NA_character_ else s
}, "")
out <- data.frame(id = smp$id, feature_a = smp$feature_a,
                  feature_b = smp$feature_b, stratum_status = smp$status,
                  live_status = live,
                  record_verdict = rec$verdict[match(smp$id, rec$id)],
                  stringsAsFactors = FALSE)
utils::write.table(out, args[2], sep = "\t", row.names = FALSE,
                   quote = TRUE, qmethod = "double")
cat("pairs:", nrow(out), "; no live row:", sum(is.na(out$live_status)),
    "; live differs from stratum:",
    sum(out$live_status != out$stratum_status, na.rm = TRUE), "\n")
d <- out[is.na(out$live_status) | out$live_status != out$stratum_status, ]
for (i in seq_len(nrow(d))) {
  cat(sprintf("%s %-22s x %-16s %-11s -> %-11s record: %s\n", d$id[i],
              d$feature_a[i], d$feature_b[i], d$stratum_status[i],
              d$live_status[i], d$record_verdict[i]))
}
