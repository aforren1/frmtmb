# Reviewer: dump brms 2.23.0 source of the functions the lane ports or
# depends on, to dev/postfit2-rev-brmssrc/<name>.R, read rather than
# recalled.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("brms")
cat("brms", format(packageVersion("brms")), "\n")
out <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2/dev/postfit2-rev-brmssrc"
dir.create(out, showWarnings = FALSE)
nms <- c("conditional_effects.brmsfit", "prepare_conditions",
         "conditional_effects.brmsterms", "conditional_effects.btl",
         "conditional_effects.mvbrmsterms", "make_point_frame",
         "validate_newdata", "prepare_predictions.brmsfit",
         "prepare_predictions_re", "prepare_predictions_re_global",
         "get_new_rsamples", "prepare_Z", "validate_re_formula",
         "conditional_smooths.brmsfit", "conditional_smooths.btl",
         "conditional_smooths.brmsterms", "make_conditions",
         "update_adterms", "posterior_average.brmsfit", "rows2labels",
         "get_group_vars", "prepare_re_newdata", "update_re_terms",
         "get_int_vars", "is_like_factor", "round_largest_remainder",
         "model_weights.brmsfit", "prepare_predictions_ranef",
         "validate_newdata_ranef", "get_levels", "new_levels_msg")
for (n in nms) {
  if (!exists(n, envir = ns, inherits = FALSE)) {
    cat("absent:", n, "\n")
    next
  }
  writeLines(deparse(get(n, envir = ns)), file.path(out, paste0(n, ".R")))
}
cat("ls matching new level / NA group:\n")
print(grep("new_lev|newlev|new_rsamp|re_newdata|rsample|na_group|ranef", ls(ns, all.names = TRUE), value = TRUE))
