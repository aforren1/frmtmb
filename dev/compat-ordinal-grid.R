# Cross every core addition term with the ordinal-only structures and
# print each pair's status and the head of its note, from a compat.tsv
# written by dev/equiv-dump-compat.R.
#   Rscript dev/compat-ordinal-grid.R <compat.tsv>
args <- commandArgs(trailingOnly = TRUE)
cp <- utils::read.delim(args[1], stringsAsFactors = FALSE)
cp <- cp[cp$in_core_only, ]
ft <- utils::read.delim(file.path(dirname(args[1]), "features.tsv"),
                        stringsAsFactors = FALSE)
aterms <- ft$name[ft$kind == "aterm"]
ord <- c("disc", "equidistant", "sum_to_zero", "cs_pred()", "thres()")
get_row <- function(a, b) {
  cp[(cp$feature_a == a & cp$feature_b == b) |
       (cp$feature_a == b & cp$feature_b == a), , drop = FALSE]
}
for (a in aterms) for (b in ord) {
  if (a == b) next
  h <- get_row(a, b)
  if (nrow(h) == 0L) next
  cat(sprintf("%-14s x %-12s %-11s %s\n", a, b, h$status[1],
              substr(h$note[1], 1, 90)))
}
