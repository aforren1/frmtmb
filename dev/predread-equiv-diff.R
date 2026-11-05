# Verdict changes of the pairs this lane remeasured, against the baseline of
# 0.78.0 (rellib-r16). Prints the changes and the two headline rates.
#   Rscript dev/predread-equiv-diff.R > dev/predread-log/equiv-diff.txt
base <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out"
out <- "dev/predread-equiv-out"
a <- read.delim(file.path(base, "results.tsv"), stringsAsFactors = FALSE)
b <- read.delim(file.path(out, "results.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(a$id, b$id))
reached <- a$id
cat("pairs:", nrow(b), " remeasured:", length(reached), "\n")
ch <- a$verdict != b$verdict | a$level %in% NA != b$level %in% NA |
  (!is.na(a$level) & !is.na(b$level) & a$level != b$level) |
  (!is.na(a$agree) & !is.na(b$agree) & a$agree != b$agree)
cat("verdict or level changes:", sum(ch), "\n")
print(data.frame(id = a$id[ch], a = a$feature_a[ch], b = a$feature_b[ch],
                 was = a$verdict[ch], now = b$verdict[ch],
                 lvl_was = a$level[ch], lvl_now = b$level[ch]),
      row.names = FALSE)
# the metric of each reached pair at L2, was and now
r <- a$id %in% reached
m <- data.frame(id = a$id[r], pair = paste(a$feature_a[r], "x", a$feature_b[r]),
                verdict = b$verdict[r], l2_was = a$l2_metric[r],
                l2_now = b$l2_metric[r], l2_what = substr(b$l2_what[r], 1, 40))
print(m[!is.na(m$l2_now) | !is.na(m$l2_was), ], row.names = FALSE)
for (d in c(base, out)) {
  s <- readLines(file.path(d, "summary.md"))
  cat(d, ":\n", s[grepl("agree with brms [(]the headline", s)], "\n")
}
vm <- a$verdict_measured != b$verdict_measured
cat("measured verdict changes (before the adjudication table):", sum(vm), "\n")
print(data.frame(id = a$id[vm], pair = paste(a$feature_a[vm], "x", a$feature_b[vm]),
                 was = a$verdict_measured[vm], now = b$verdict_measured[vm]),
      row.names = FALSE)
