# Fix round 1 of lane catmulti: the new tests seen to fail on the
# lane's ROUND-1 behavior. Each fix is undone in this session only, by
# putting its round-1 body back into the loaded namespace (no library
# changes), and the test files are run against that.
#   Rscript dev/catmulti-revert-arm.R
# Output: dev/catmulti-out/revert-arm.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(testthat); library(frmtmb); library(frmtmb.sample)
})
put <- function(pkg, nm, f) {
  ns <- asNamespace(pkg)
  environment(f) <- ns
  unlockBinding(nm, ns)
  assign(nm, f, envir = ns)
  lockBinding(nm, ns)
}
edit_body <- function(pkg, nm, old, new) {
  f <- get(nm, asNamespace(pkg))
  src <- paste(deparse(f), collapse = "\n")
  stopifnot(grepl(old, src, fixed = TRUE))
  put(pkg, nm, eval(parse(text = sub(old, new, src, fixed = TRUE))))
}
# 1. the raw response name in the lncor prefix
edit_body("frmtmb", "lncor_blocks", "brms_stan_name(r)", "r")
edit_body("frmtmb.sample", "draws_lncor_cols", "brms_stan_name(r)", "r")
put("frmtmb", "resp_from_stan_name", function(spec, resp) resp)
# 2. the df counted every outer parameter
put("frmtmb", "loglik_df", function(object) n_outer_est(object))
# 3. transform = TRUE refused only with two or more predictors
edit_body("frmtmb.sample", "posterior_linpred.frmtmb_draws",
          "if (length(stack) && isTRUE(transform))",
          "if (length(stack) > 1L && isTRUE(transform))")
# 4. no refcat = NA sentence in the SE warning
edit_body("frmtmb", "se_lost_message",
          "is.function(r$family[[\"refcat_first\"]])", "FALSE")
# 5. the categorical display refused its Wald band
put("frmtmb", "ce_fitted_bands",
    function(fam) ce_matrix_mean(fam))
files <- list(
  c("frmtmb", "tests/testthat/test-refcat.R"),
  c("frmtmb", "tests/testthat/test-simplex-mvbf.R"),
  c("frmtmb.sample",
    "extensions/frmtmb.sample/tests/testthat/test-linpred-stack.R"),
  c("frmtmb.sample",
    "extensions/frmtmb.sample/tests/testthat/test-simplex-mvbf-draws.R"))
for (pf in files) {
  r <- as.data.frame(test_file(pf[2], package = pf[1],
                               env = test_env(pf[1]), reporter = "silent",
                               stop_on_failure = FALSE))
  cat(sprintf("REVERTED %s pass=%d fail=%d err=%d\n", basename(pf[2]),
              sum(r$passed), sum(r$failed), sum(r$error)))
  bad <- r[r$failed > 0 | r$error, "test"]
  for (b in bad) cat("  not passing:", b, "\n")
}
