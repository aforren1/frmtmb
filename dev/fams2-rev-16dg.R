# Reviewer, claim 3: rows 16d to 16g of the gated brms likelihood tier,
# block by block so each LPCHECK line has its label. The four programs
# were compiled fresh into the reviewer's cache copy on the first gated
# run (the worker's four entries were removed from the copy).
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
sp <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
             "c--Users-adf44-source-r-frmtmb/",
             "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad")
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(sp, "stan-cache"))
options(frmtmb.brms_lp_report = TRUE)
suppressPackageStartupMessages(library(testthat))
f <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/tests/testthat/test-brms-likelihood.R"
for (desc in c("row 16d: extended-support beta with kappa ~ x",
               "row 16e: zero-inflated beta-binomial with zi ~ x",
               "row 16f: hurdle cumulative with hu ~ x",
               "row 16g: hurdle cumulative, probit, with disc ~ 0 + x")) {
  cat("\n####", desc, "\n")
  r <- as.data.frame(test_file(f, package = "frmtmb",
                               env = test_env("frmtmb"), reporter = "summary",
                               desc = desc, load_package = "installed"))
  cat(sprintf("RESULT [%s] pass=%d fail=%d err=%d skip=%d\n", desc,
              sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
}
