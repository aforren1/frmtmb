# Every GitHub workflow of the tree parses as YAML, and every
# .github/scripts/*.R file parses as R.
#   Rscript dev/rel077-workflows.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel077")
fs <- list.files(".github/workflows", "[.]ya?ml$", full.names = TRUE)
bad <- 0
for (f in fs) {
  r <- tryCatch({ y <- yaml::read_yaml(f); if (is.null(y$jobs)) "no jobs" else NULL },
                error = function(e) conditionMessage(e))
  if (!is.null(r)) { bad <- bad + 1; cat("BAD", f, r, "\n") }
}
cat("workflows:", length(fs), "parse errors:", bad, "\n")
ss <- list.files(".github/scripts", "[.]R$", full.names = TRUE)
pb <- sum(vapply(ss, function(f) inherits(try(parse(f), silent = TRUE), "try-error"), NA))
cat(".github/scripts R files:", length(ss), "parse errors:", pb, "\n")
