# The sampled pairs this lane's change reaches, and their baseline verdicts:
# every pair with a method reader, a cs() term, or a me() term. Copies the
# baseline of the others into dev/predread-equiv-out/, so that equiv-run.sh
# measures only these.
#   Rscript dev/predread-equiv-select.R
base <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out"
out <- "dev/predread-equiv-out"
r <- read.delim(file.path(base, "results.tsv"), stringsAsFactors = FALSE)
reach <- c("frm_curve", "frm_curve_contrast", "fitted", "predict", "simulate",
           "residuals", "emmeans", "cs_pred()", "me()", "frm_lp_basis",
           "frm_linpred", "conditional_effects")
hit <- r$feature_a %in% reach | r$feature_b %in% reach
cat("pairs:", nrow(r), " reached:", sum(hit), "\n")
print(r[hit, c("id", "feature_a", "feature_b", "verdict", "level")],
      row.names = FALSE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (id in r$id[!hit]) {
  file.copy(file.path(base, paste0(id, ".tsv")), out, overwrite = TRUE)
  rds <- file.path(base, paste0(id, ".rds"))
  if (file.exists(rds)) file.copy(rds, out, overwrite = TRUE)
}
writeLines(r$id[hit], "dev/predread-log/equiv-reached.txt")
