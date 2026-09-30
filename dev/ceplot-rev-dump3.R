# Reviewer helper (lane ceplot, punch round 1): deparse the functions
# whose mutant patterns changed.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
for (f in c("ce_valid_effects", "ce_model_vars", "new_level_pick_apply",
            "old_level_pick_add", "ps_select")) {
  ns <- if (f == "ps_select") "frmtmb.sample" else "frmtmb"
  cat("#####", f, "\n")
  writeLines(deparse(get(f, asNamespace(ns))))
}
