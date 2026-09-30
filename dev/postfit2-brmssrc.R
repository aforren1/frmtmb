# Dumps brms 2.23.0's source for the functions this lane ports, from
# the installed namespace, into dev/postfit2-brmssrc/.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
stopifnot(packageVersion("brms") == "2.23.0")
ns <- asNamespace("brms")
out <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2/dev/postfit2-brmssrc"
dir.create(out, showWarnings = FALSE)
pat <- c("conditional_smooths", "conditional_effects", "make_conditions",
         "posterior_average", "update_adterms", "spaghetti", "too_far",
         "select_points", "make_point_frame", "prepare_conditions",
         "conditional_effects_internal", "make_surface",
         "get_int_vars", "rows2labels", "validate_conditions",
         "get_cond__", "stanfit_", "average", "model_weights",
         "get_all_effects", "sm_labels", "get_se_", "ad_terms",
         "allow_adterms", "fix_", "posterior_samples")
nms <- ls(ns, all.names = TRUE)
sel <- unique(unlist(lapply(pat, function(p) grep(p, nms, value = TRUE))))
for (nm in sel) {
  f <- get(nm, envir = ns)
  if (!is.function(f)) next
  writeLines(c(paste0(nm, " <- "), deparse(f, control = "useSource")),
             file.path(out, paste0(gsub("[^A-Za-z0-9_.]", "_", nm), ".R")))
}
cat(sort(sel), sep = "\n")
