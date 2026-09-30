# Reviewer helper (lane ceplot): deparse the functions the mutants of
# dev/ceplot-rev-mutants.R edit, so their patterns match the installed
# code exactly.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
for (f in c("ce_plan_eval", "ce_valid_effects", "new_level_pick_apply",
            "old_level_pick_add", "ce_mm_in_block")) {
  cat("#####", f, "\n")
  writeLines(deparse(get(f, asNamespace("frmtmb"))))
}
cat("##### frmtmb.sample imports of core internals\n")
imp <- parent.env(asNamespace("frmtmb.sample"))
print(intersect(ls(imp), c("ce_plan_eval", "ce_level_plan", "ce_plan_part",
                           "ce_valid_effects", "ce_grids_build",
                           "new_level_pick_apply")))
