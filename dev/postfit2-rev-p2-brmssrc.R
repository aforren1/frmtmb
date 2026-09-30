# Reviewer, punch round 2: brms 2.23.0 source of how new group levels
# get their J indices (to judge the mm() shared-draw behavior).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("brms")
hits <- Filter(function(n) {
  f <- get(n, envir = ns)
  is.function(f) && any(grepl("new_levels|new levels|max_level|expand_levels", deparse(f)))
}, ls(ns, all.names = TRUE))
print(hits)
for (n in intersect(hits, c("data_gr_local", ".data_gr_local", "standata_basis_gr", "get_group_vars", "subset_levels", "prepare_predictions_re_global", "data_re", ".data_re"))) {
  cat("=====", n, "\n"); print(get(n, envir = ns))
}
