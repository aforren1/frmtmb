# Reviewer, lane sampfix, script 08: how brms 2.23.0's conditional_effects()
# treats the grouping variable at re_formula = NULL (read from source).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
for (fn in c("conditional_effects.brmsfit", "conditional_effects.brmsterms",
             "prepare_conditions")) {
  f <- deparse(get(fn, asNamespace("brms")))
  cat("==", fn, "\n")
  print(grep("new_levels|re_formula|sample_new|NA_", f, value = TRUE))
}
