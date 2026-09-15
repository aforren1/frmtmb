root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L)
    stop("many in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

# BLOCKER A: the hooks and the unload hook go with the binding swap.
sub1("R/zzz.R",
paste0("#' Load hook for the optional downstream integrations. It first hands\n",
       "#' every shared generic back to the package that owns it, which has to\n",
       "#' happen here rather than in the NAMESPACE because the owners are\n",
       "#' optional (see `R/generic-owners.R`). It then registers the\n"),
paste0("#' Load hook for the optional downstream integrations. It first turns\n",
       "#' every borrowed generic into an active binding that resolves to the\n",
       "#' package that owns it, which has to happen here rather than in the\n",
       "#' NAMESPACE because the owners are optional, and here rather than\n",
       "#' anywhere else because the namespace is still unsealed (see\n",
       "#' `R/generic-owners.R`). It then registers the\n"))

sub1("R/zzz.R",
"  frm_adopt_generics(pkgname)\n",
"  frm_install_generics(pkgname)\n")

sub1("R/zzz.R",
paste0("\n#' Unload hook. The generic-adoption load hooks point into this\n",
       "#' namespace, so they come off when it goes away.\n",
       "#'\n",
       "#' @noRd\n",
       ".onUnload <- function(libpath) {\n",
       "  frm_drop_generic_hooks()\n",
       "  invisible()\n",
       "}\n"),
"")

# N9: the comment describes the defect this change removed.
sub1("R/loo.R",
paste0("  # frmtmb's generic would otherwise mask loo's own function for anyone\n",
       "  # who attaches both, and `loo_compare(loo(d1), loo(d2))`, the\n",
       "  # spelling the sampling package's help page recommends, would stop at\n",
       "  # \"no applicable method\", because loo does not export its default\n",
       "  # method for the search path to find. Reaching for loo's METHOD rather\n",
       "  # than calling loo::loo_compare() is deliberate: dispatch from inside\n",
       "  # this namespace would find this function again and recurse forever.\n"),
paste0("  # Reached only when loo is NOT loaded, because frmtmb's exported\n",
       "  # `loo_compare` is an active binding that resolves to loo's own\n",
       "  # generic whenever loo is there (R/generic-owners.R). Before that\n",
       "  # change this method carried the whole burden: frmtmb's generic\n",
       "  # masked loo's function for anyone who attached both, and\n",
       "  # `loo_compare(loo(d1), loo(d2))` stopped at \"no applicable\n",
       "  # method\" because loo does not export its default method for the\n",
       "  # search path to find. Reaching for loo's METHOD rather than\n",
       "  # calling loo::loo_compare() is still deliberate: dispatch from\n",
       "  # inside this namespace would find this function again and recurse.\n"))
cat("DONE\n")
