# Emit the per-vignette markdown tables for dev/brms-vignette-port.md.
HERE <- "C:/Users/adf44/source/r/frmtmb-wt-audit/dev/brms-port"
source(file.path(HERE, "port-lib.R"))
M <- readRDS(file.path(HERE, "results-merged.rds"))
VIGS <- c("brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
          "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
          "brms_missings", "brms_customfamilies")
short <- function(s, n) gsub("\\|", "\\\\|", substr(gsub("\\s+", " ", s), 1, n))
cls <- classify
for (v in VIGS) {
  s <- M[M$vignette == v & M$kind %in% c("model", "post"), ]
  cat("\n#### ", v, "\n\n", sep = "")
  cat("| id | kind | code | class | note |\n|---|---|---|---|---|\n")
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    note <- if (r$status_spell == "OK") "" else short(r$msg_spell, 90)
    cat("| ", sub("^[a-z_]+\\.", "", r$id), " | ", r$kind, " | `", short(r$src, 62),
        "` | ", cls(r), " | ", note, " |\n", sep = "")
  }
}
