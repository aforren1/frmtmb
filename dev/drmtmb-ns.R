# Which of the features drmTMB's website names exist in the CRAN 0.7.0
# namespace, exported or not.
.libPaths(c("C:/Users/adf44/source/r/drmtmb-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- ls(getNamespace("drmTMB"), all.names = TRUE)
cat("drmTMB", format(packageVersion("drmTMB")), ":", length(ns),
    "objects in namespace\n")
pats <- c("herit", "icc", "repeatab", "chibar", "lrt", "boundary_p",
          "objective_at", "worm", "centile", "phylo_interaction", "coevol",
          "meta_V", "sample", "tmbstan", "mcmc")
for (pat in pats) {
  hit <- head(grep(pat, ns, value = TRUE, ignore.case = TRUE), 12)
  cat(pat, ":", paste(hit, collapse = ", "), "\n")
}
cat("\nexports:\n")
print(sort(getNamespaceExports("drmTMB")))
cat("\nanova.drmTMB:\n")
print(drmTMB:::anova.drmTMB)
uses_pchisq <- Filter(function(nm) {
  f <- get(nm, getNamespace("drmTMB"))
  is.function(f) && any(grepl("pchisq", deparse(f), fixed = TRUE))
}, ns)
cat("functions calling pchisq:", length(uses_pchisq), "\n")
