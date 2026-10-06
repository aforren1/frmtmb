# Like dev/nanse-run-tests.R, but every firing of the standard-error
# warning also records what the base build would have reported for the
# same fit (sdreport()'s own cov.fixed: how many SEs are non-finite,
# whether pdHess) and the spectrum the verdict was read from. The
# objective's state is saved and restored around the extra sdreport().
#   Rscript dev/nanse-run-diag.R <package> <test file> <fire log>
LIB <- "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
fire <- a[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", find.package("frmtmb"), " | ", find.package(p), "\n", sep = "")
ns <- asNamespace("frmtmb")
.nanse_diag <- function(fit, msg, fire) {
  ns <- asNamespace("frmtmb")
  st <- ns$obj_state_save(fit$obj)
  on.exit(ns$obj_state_restore(fit$obj, st))
  base <- tryCatch({
    s <- suppressWarnings(RTMB::sdreport(fit$obj))
    se <- suppressWarnings(sqrt(diag(s$cov.fixed)))
    sprintf("base: %d of %d SE non-finite, pdHess %s", sum(!is.finite(se)),
            length(se), s$pdHess)
  }, error = function(e) paste("base sdreport error:", conditionMessage(e)))
  an <- fit$cache$hessian_fixed$an
  H <- fit$cache$hessian_fixed$H
  ev <- if (!is.null(an$ev)) signif(an$ev, 3) else NA
  re <- length(fit$obj$env$random) > 0
  cat(gsub("[\r\n\t]+", " ", msg), "\n  call: ",
      substr(deparse1(fit$call), 1, 200), "\n  ", base,
      " | random effects ", re,
      " | code ", fit$opt$convergence, " | maxgrad ",
      signif(max(abs(fit$obj$gr(fit$opt$par))), 3),
      "\n  scaled eigenvalues: ", paste(ev, collapse = " "),
      "\n  diag H: ", if (is.null(H)) "NA" else {
        paste(signif(diag(H), 3), collapse = " ")
      },
      # the quantity se_empty_row is compared with, per parameter, and
      # whether the rule emptied that row
      "\n  rowmax: ", if (is.null(H)) "NA" else {
        u <- fit$par_units %||% rep(1, nrow(H))
        rm <- apply(abs(H * outer(u, u)), 1L, max)
        lost <- fit$cache$hessian_fixed$an$lost
        nm <- ns$outer_par_names(fit)
        paste(sprintf("%s=%.3g%s", nm, rm,
                      ifelse(nm %in% names(lost), "*", "")), collapse = " ")
      }, "\n",
      file = fire, append = TRUE, sep = "")
}
assign(".nanse_diag", .nanse_diag, envir = globalenv())
assign(".nanse_fire", fire, envir = globalenv())
suppressMessages(trace(
  "se_lost_message", where = ns, print = FALSE,
  exit = quote(get(".nanse_diag", envir = globalenv())(
    fit, returnValue(), get(".nanse_fire", envir = globalenv())))))
r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
