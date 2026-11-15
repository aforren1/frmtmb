# Item 5516: split the recorded dev/equiv-data/results.tsv into the
# per-pair files dev/equiv-summary.R reads, so the summary can be rerun
# into a scratch outdir (never --record) to show the new count line.
#   Rscript dev/hygiene-equiv-split.R dev/hygiene-out/equiv-sum
out <- commandArgs(trailingOnly = TRUE)[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
r <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE,
                na.strings = "NA")
for (i in seq_len(nrow(r))) {
  utils::write.table(r[i, ], file.path(out, paste0(r$id[i], ".tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}
cat("pairs", nrow(r), "\n")
k <- r$verdict == "refusal correct" & r$status == "conditional"
cat("from the record directly: refusal correct on a conditional row",
    sum(k), "unique,", sum(r$draws[k]), "draws\n")
