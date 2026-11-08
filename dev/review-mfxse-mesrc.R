# Review of lane mfxse: dump marginaleffects 1.0.0's functions that call
# get_coef(), get_vcov(), set_coef(), or read hypothesis_is_formula, to
# check the signal the lane keys hypotheses() off.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
cat("marginaleffects", format(packageVersion("marginaleffects")), "\n")
out <- "C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/review-mfxse-out/mesrc"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns)
  if (!is.function(f)) next
  txt <- paste(deparse(f, control = "useSource"), collapse = "\n")
  hit <- grepl("hypothesis_is_formula|get_coef\\(|get_vcov\\(|set_coef\\(|get_jacobian|numDeriv|numderiv", txt)
  if (hit) {
    writeLines(txt, file.path(out, paste0(gsub("[^A-Za-z0-9_.]", "_", nm), ".R")))
    cat(nm, "\n")
  }
}
