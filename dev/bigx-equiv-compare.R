# Lane bigx: each pair's measured verdict, level and agreement on the
# lane build (dev/bigx-equiv-out) against the survey's record
# (dev/equiv-data/results.tsv, measured on 0.79.0). The record's
# `verdict` is adjudicated; `verdict_measured` is what the harness
# itself says, and that is what a remeasure gives.
#   Rscript dev/bigx-equiv-compare.R
rec <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
fs <- list.files(Sys.getenv("BIGX_EQUIV_OUT", "dev/bigx-equiv-out"), "[.]tsv$", full.names = TRUE)
new <- do.call(rbind, lapply(fs, function(f) {
  read.delim(f, stringsAsFactors = FALSE)[, c("id", "verdict", "level",
                                               "agree", "l1_const",
                                               "frmtmb")]
}))
cat("pairs in the record", nrow(rec), "; remeasured", nrow(new), "\n")
cat("lane build version read:", unique(new$frmtmb), "\n")
m <- merge(rec[, c("id", "verdict_measured", "verdict", "level", "agree",
                   "l1_const")],
           new, by = "id", suffixes = c(".rec", ".lane"))
neq <- function(a, b) !(is.na(a) & is.na(b)) & (is.na(a) | is.na(b) | a != b)
ch <- m[which(neq(m$verdict_measured, m$verdict.lane) |
                neq(m$level.rec, m$level.lane) |
                neq(m$agree.rec, m$agree.lane)), ]
cat("pairs whose measured verdict, level or agreement changed:", nrow(ch),
    "\n")
if (nrow(ch)) {
  print(ch[, c("id", "verdict_measured", "verdict.rec", "verdict.lane",
               "level.rec", "level.lane", "agree.rec", "agree.lane",
               "l1_const.rec", "l1_const.lane")], row.names = FALSE)
}
missing <- setdiff(rec$id, new$id)
cat("not remeasured:", if (length(missing)) missing else "none", "\n")
