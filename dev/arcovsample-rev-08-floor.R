# REVIEW script 08, claim 7: the worker's frmtmb.sample against the
# REFERENCE frmtmb. The floor must fail LOUDLY.
#
#   Rscript dev/arcovsample-rev-08-floor.R
#
# frmtmb.sample's NAMESPACE has import(frmtmb), not an explicit
# importFrom for the two new names, so nothing stops the install or the
# load; the question is what happens at the call.

FLOOR <- "C:/Users/adf44/source/r/wt-arcovsample-rev-floorlib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(FLOOR, REF, USER))
Sys.setenv(NOT_CRAN = "true")
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}

cat("--- does the package even load?\n")
lr <- tryCatch({ suppressMessages(library(frmtmb.sample)); "loaded" },
               error = function(e) paste("LOAD ERROR:",
                                         conditionMessage(e)))
cat("  ", lr, "\n", sep = "")
cat("  frmtmb        ", format(packageVersion("frmtmb")), " at ",
    dirname(system.file("DESCRIPTION", package = "frmtmb")), "\n", sep = "")
cat("  frmtmb.sample ", format(packageVersion("frmtmb.sample")), " at ",
    dirname(system.file("DESCRIPTION", package = "frmtmb.sample")), "\n",
    sep = "")
cat("  reference frmtmb exports arma_cond_resp: ",
    "arma_cond_resp" %in% getNamespaceExports("frmtmb"), "\n", sep = "")

cat("\n--- log_lik() on a cov = FALSE ARMA fit\n")
set.seed(13L)
dd <- data.frame(g = factor(rep(1:4, each = 6L)), t = rep(1:6, 4L))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(nrow(dd), 0, 0.8)
fit <- frm(bf(y ~ x + ar(t, g)), family = gaussian(), data = dd,
           dry_run = "objective")
lab <- c(frmtmb::brms_par_labels(fit), "lp__")
m <- matrix(0, 3L, length(lab), dimnames = list(NULL, lab))
m[, "b_Intercept"] <- 0.5; m[, "b_x"] <- 0.4
m[, "sigma"] <- 0.8; m[, "thetaac_1"] <- 0.4
ds <- structure(list(stanfit = NULL, draws = m, fit = fit),
                class = "frmtmb_draws")
r <- tryCatch(log_lik(ds), error = identity, warning = identity)
if (inherits(r, "condition")) {
  cat("  ", class(r)[1L], ": ", conditionMessage(r), "\n", sep = "")
} else {
  cat("  NO ERROR. dim = ", paste(dim(r), collapse = " x "),
      "  any non-finite: ", any(!is.finite(r)), "\n", sep = "")
  cat("  SILENT: rowSums = ", paste(format(rowSums(r), digits = 8),
                                    collapse = " "), "\n", sep = "")
}

cat("\n--- posterior_predict() on the same fit (the other new call)\n")
r2 <- tryCatch({ set.seed(1); posterior_predict(ds) }, error = identity)
if (inherits(r2, "condition")) {
  cat("  ", class(r2)[1L], ": ", conditionMessage(r2), "\n", sep = "")
} else {
  cat("  NO ERROR, dim ", paste(dim(r2), collapse = " x "), "\n", sep = "")
}

cat("\n--- the whole test file\n")
suppressMessages(library(testthat))
res <- tryCatch(
  as.data.frame(test_file(
    "extensions/frmtmb.sample/tests/testthat/test-loo.R",
    package = "frmtmb.sample",
    env = testthat::test_env("frmtmb.sample"), reporter = "silent")),
  error = function(e) { cat("  LOADERROR ", conditionMessage(e), "\n",
                            sep = ""); NULL })
if (!is.null(res)) {
  cat("RESULT floor test-loo.R pass=", sum(res$passed), " fail=",
      sum(res$failed), " err=", sum(res$error), " skip=",
      sum(res$skipped), "\n", sep = "")
  bad <- res[res$failed > 0 | res$error > 0, , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    cat("--- BAD: ", bad$test[i], "\n", sep = "")
    for (rs in bad$result[[i]]) {
      if (inherits(rs, c("expectation_failure", "expectation_error"))) {
        cat("    ", conditionMessage(rs), "\n", sep = "")
      }
    }
  }
}
cat("\nDONE\n")
