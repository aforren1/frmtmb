# Reviewer: list every FAIL/CASCADE row of the mechanical spell pass
# with its spell-pass message, so the lane's root-cause charging can be
# recounted by hand for chosen causes.
#
#   Rscript dev/vigport-rev-gaps.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
M <- readRDS(file.path(root, "vigport-port-out/r5/results-merged.rds"))
bad <- M[M$bucket %in% c("FAIL", "CASCADE"), ]
cat("rows:", nrow(bad), " by kind:",
    paste(names(table(bad$kind)), table(bad$kind), collapse = ", "), "\n\n")
for (i in seq_len(nrow(bad))) {
  cat(sprintf("%-26s %-5s %-8s %s\n      src: %s\n", bad$id[i], bad$kind[i],
              bad$bucket[i], substr(bad$msg_spell[i], 1, 90),
              substr(bad$src_run[i], 1, 90)))
}
