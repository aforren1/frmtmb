# Review: one core test file with every extension loaded first.
#   Rscript dev/review-rel075-preload1.R <file>
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(testthat); library(frmtmb)})
for (q in c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
            "frmtmb.ode", "frmtmb.sample", "frmtmb.survival", "frmtmb.spatial"))
  loadNamespace(q)
cat("loaded:", paste(sort(loadedNamespaces()[startsWith(loadedNamespaces(), "frmtmb")]), collapse = " "), "\n")
r <- as.data.frame(test_file(commandArgs(TRUE)[1], package = "frmtmb",
  env = testthat::test_env("frmtmb"), reporter = "summary"))
cat("RESULT pass=", sum(r$passed), " fail=", sum(r$failed), " err=", sum(r$error),
    " skip=", sum(r$skipped), " warn=", sum(r$warning), "\n", sep = "")
