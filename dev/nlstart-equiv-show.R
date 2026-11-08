# Print level/agree/reason of pairs in a remeasure dir beside the record.
# Usage: Rscript dev/nlstart-equiv-show.R <dir> <id> ...
a <- commandArgs(TRUE)
rec <- utils::read.delim("dev/equiv-data/results.tsv", colClasses = "character")
for (id in a[-1]) {
  now <- utils::read.delim(file.path(a[1], paste0(id, ".tsv")),
                           colClasses = "character")
  r <- rec[rec$id == id, ]
  for (k in c("level", "agree", "reason", "frm_ok", "l1_const", "frm_warn")) {
    cat(id, k, "| record:", substr(r[[k]], 1, 150), "| now:",
        substr(now[[k]], 1, 150), "\n")
  }
}
