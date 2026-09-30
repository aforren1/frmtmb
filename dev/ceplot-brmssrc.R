# Dumps brms 2.23.0's source for the functions lane ceplot ports, from
# the installed namespace, into dev/ceplot-brmssrc/.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
stopifnot(packageVersion("brms") == "2.23.0")
ns <- asNamespace("brms")
out <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-brmssrc"
dir.create(out, showWarnings = FALSE)
pat <- c("^plot[.]", "conditional_effects", "^fitted", "posterior_samples",
         "^parnames", "^nsamples", "^nuts_params", "^hypothesis",
         "^plot_hyp", "^get_all_effects", "prepare_conditions",
         "make_point_frame", "^\\.plot", "^ndraws", "^variables",
         "^posterior_epred", "prepare_predictions", "old_levels",
         "new_levels", "^get_dpar", "^conditional_smooths", "^mi_",
         "^validate_newdata", "^rug", "^theme")
nms <- ls(ns, all.names = TRUE)
sel <- unique(unlist(lapply(pat, function(p) grep(p, nms, value = TRUE))))
for (nm in sel) {
  f <- get(nm, envir = ns)
  if (!is.function(f)) next
  writeLines(c(paste0(nm, " <- "), deparse(f, control = "useSource")),
             file.path(out, paste0(gsub("[^A-Za-z0-9_.]", "_", nm), ".R")))
}
cat(sort(sel), sep = "\n")
