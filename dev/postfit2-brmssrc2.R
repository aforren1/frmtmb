# Second dump of brms 2.23.0 helpers the ported functions call.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
stopifnot(packageVersion("brms") == "2.23.0")
ns <- asNamespace("brms")
out <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2/dev/postfit2-brmssrc"
sel <- c("validate_weights", "validate_weights_method", "split_dots",
         "round_largest_remainder", "posterior_smooths",
         "posterior_smooths.brmsfit", "posterior_smooths.btl",
         "fill_newdata", "add_effects__", "prepare_cond_data",
         "scale_unit", "is_like_factor", "frame_sm", "exclude_terms",
         "exclude_terms.brmsfit", "exclude_terms.brmsformula",
         "exclude_terms.mvbrmsformula", "combine_prefix",
         "posterior_summary", "posterior_summary.default",
         "validate_draw_ids", "as_formula", "formula2str", "get_matches",
         "deparse0", "loo_model_weights.brmsfit", "pp_average.brmsfit")
for (nm in sel) {
  f <- get(nm, envir = ns)
  writeLines(c(paste0(nm, " <- "), deparse(f, control = "useSource")),
             file.path(out, paste0(nm, ".R")))
}
