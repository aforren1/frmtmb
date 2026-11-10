# Reviewer scratch (review of lane translator, 2026-11-09): the control
# of review-translator-p048fix.R for p047, at the SURVEY points (no
# l1_points). Was: the remedy
# for the blind t2() blocks of p048 and p047. Read L1 at the lane's own
# fixed points (eq_jitter_optimizer(), frm()'s start moved by 0.3 sin),
# where every smooth SD sits near its start rather than at ~1e-4, and
# rerun the one-block map mutant that the survey's points could not see.
# Records nothing.
#
#   Rscript dev/review-translator-p048fix.R
lib <- "C:/Users/adf44/source/r/rellib-r19"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             winslash = "/"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
  requireNamespace("brms")
  requireNamespace("rstan")
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
fmt <- function(x) if (!length(x) || !is.numeric(x)) "-" else
  paste(signif(x, 3), collapse = ",")
fresh <- function() {
  e <- new.env(parent = asNamespace("frmtmb"))
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
    sys.source(f, envir = e)
  }
  e
}
smp <- read.delim("dev/equiv-data/sample.tsv")
at_points <- function(s) {
  return(s)
  s$l1_points <- list(scale = 0.3, why = "review: SDs away from 0")
  s
}
run <- function(label, id, block = NULL) {
  e <- fresh()
  row <- smp[smp$id == id, ]
  key <- paste(row$feature_a, row$feature_b, sep = " x ")
  mk <- e$eq_specs[[key]]
  e$eq_specs[[key]] <- function() at_points(mk())
  if (!is.null(block)) {
    v <- e$eq_rules$zs_prefixed$value
    e$eq_rules$zs_prefixed$value <- function(nm, ...) {
      out <- v(nm, ...)
      if (identical(nm, block)) {
        out$value <- array(rev(out$value), length(out$value))
      }
      out
    }
  }
  r <- e$eq_measure(key, row)
  cat(sprintf("%-44s %s verdict %-12s L1 ok %-5s const %s grad %s err %s\n",
              label, id, r$verdict %||% "", format(r$l1_ok %||% NA),
              fmt(c(r$l1_const, r$l1_extra_const)), fmt(r$l1_grad),
              substr(r$l1_error %||% "", 1, 100)))
}




run("survey points: zs_b1_1_3 alone reversed", "p047", "zs_b1_1_3")
