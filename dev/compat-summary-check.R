# Rebuild per-pair tsvs from the record (dev/equiv-data/results.tsv,
# measured verdicts, summary columns dropped) into <outdir>, then overlay
# the pairs measured by the edited dev/equiv-measure.R from the given
# directories. dev/equiv-summary.R then runs on <outdir> without --record.
#   Rscript dev/compat-summary-check.R <outdir> [--live=<tsv>] [overlay ...]
# --live gives the rebuilt rows the live status of dev/equiv-live-row.R's
# output, as dev/equiv-measure.R would record it (an overlaid pair keeps
# its own).
args <- commandArgs(trailingOnly = TRUE)
live_arg <- grep("^--live=", args, value = TRUE)
args <- setdiff(args, live_arg)
out <- args[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rec <- utils::read.delim("dev/equiv-data/results.tsv",
                         stringsAsFactors = FALSE, colClasses = "character")
# the per-pair columns are those eq_flatten() writes; the summary adds
# the rest
drop <- c("verdict_measured", "adjudicated", "top_level", "weak_evidence",
          "measured", "nodefect", "equivalent")
rec$verdict <- rec$verdict_measured
rec <- rec[setdiff(names(rec), drop)]
if (length(live_arg)) {
  lv <- utils::read.delim(sub("^--live=", "", live_arg),
                          stringsAsFactors = FALSE)
  rec$stratum_status <- rec$status
  rec$live_status <- lv$live_status[match(rec$id, lv$id)]
}
for (i in seq_len(nrow(rec))) {
  utils::write.table(rec[i, ], file.path(out, paste0(rec$id[i], ".tsv")),
                     sep = "\t", row.names = FALSE, quote = TRUE,
                     qmethod = "double")
}
n_ov <- 0L
for (d in args[-1]) {
  for (f in list.files(d, "[.]tsv$", full.names = TRUE)) {
    file.copy(f, file.path(out, basename(f)), overwrite = TRUE)
    n_ov <- n_ov + 1L
  }
}
cat("rebuilt", nrow(rec), "pairs; overlaid", n_ov, "\n")
