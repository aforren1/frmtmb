# Parse every workflow and site YAML file of the release tree, and print
# the values the merge touched: the steps of the merged pkgdown job, the
# gate steps of R-CMD-check.yaml and test-coverage.yaml, and the files
# the helper identity job hashes.
#   Rscript dev/rel075-yaml.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel075")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
fs <- c(list.files(".github/workflows", "[.]ya?ml$", full.names = TRUE),
        "_pkgdown.yml",
        list.files("extensions", "^_pkgdown[.]yml$", recursive = TRUE,
                   full.names = TRUE))
bad <- 0
for (f in fs) {
  r <- tryCatch({ yaml::read_yaml(f); NULL },
                error = function(e) conditionMessage(e))
  if (!is.null(r)) { bad <- bad + 1; cat("YAML ERROR", f, r, "\n") }
}
cat("parsed", length(fs), "YAML files, errors", bad, "\n")
pk <- yaml::read_yaml(".github/workflows/pkgdown.yaml")
for (j in names(pk$jobs)) {
  cat("pkgdown job", j, "steps:\n")
  for (s in pk$jobs[[j]]$steps) cat("  -", s$name %||% s$uses, "\n")
}
cat("pkgdown push paths:", paste(pk$on$push$paths, collapse = " "), "\n")
for (w in c("R-CMD-check.yaml", "test-coverage.yaml")) {
  y <- yaml::read_yaml(file.path(".github/workflows", w))
  for (j in names(y$jobs)) {
    g <- Filter(function(s) grepl("skip-gate", s$run %||% ""), y$jobs[[j]]$steps)
    cat(w, j, "gate steps:", length(g),
        if (length(g)) paste("|", gsub("\n", " ", g[[1]]$run)), "\n")
  }
}
