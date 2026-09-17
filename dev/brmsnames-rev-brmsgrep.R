## Reviewer: print brms internals that decide names.
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("brms")
for (f in c("rename", "change_re", "change_re_levels", "get_cornames",
            "combine_prefix", "rename_pars", "make_index_names",
            "hypothesis_internal", "eval_hypothesis", "replace_special_chars",
            "validate_resp", "get_rnames", "parse_resp", "prepare_rename")) {
  if (exists(f, ns)) { cat("=====", f, "\n"); print(get(f, ns)) }
}
hits <- character(0)
for (n in ls(ns)) {
  o <- get(n, ns)
  if (is.function(o)) {
    d <- paste(deparse(o), collapse = "\n")
    if (grepl("levels", d) && grepl('"_"', d) && grepl("gsub", d))
      hits <- c(hits, n)
  }
}
cat("functions with gsub, levels and underscore:", hits, "\n")
