# Reviewer scratch (second fix-round check of lane translator): spot
# check of two pairs the per-row audit called blind because the
# correlation sits at +/-1 (diag(L_1)[2] 0.0017 for p053, 0.0038 for
# p116). Row 2 of z_1 alone is reversed after the translator, at the
# survey's points (the spec's l1_points removed) and at the lane's fixed
# points (eq_degenerate). Records nothing.
#
#   Rscript dev/review-translator-corrow.R
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
smp <- read.delim("dev/equiv-data/sample.tsv")
run <- function(id, label, rev_row = NULL, survey = FALSE) {
  row <- smp[smp$id == id, ]
  key <- paste(row$feature_a, row$feature_b, sep = " x ")
  e <- new.env(parent = asNamespace("frmtmb"))
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
    sys.source(f, envir = e)
  }
  if (survey) {
    mk <- e$eq_specs[[key]]
    e$eq_specs[[key]] <- function() {
      s <- mk()
      s$l1_points <- NULL
      s
    }
  }
  if (!is.null(rev_row)) {
    tr <- e$stan_pars_from_fit
    e$stan_pars_from_fit <- function(...) {
      p <- tr(...)
      z <- p[["z_1"]]
      z[rev_row, ] <- rev(z[rev_row, ])
      p[["z_1"]] <- z
      p
    }
  }
  r <- e$eq_measure(key, row)
  cat(sprintf("%s %-36s verdict %-12s L1 ok %-5s const %s grad %s err %s\n",
              id, label, r$verdict %||% "", format(r$l1_ok %||% NA),
              fmt(c(r$l1_const, r$l1_extra_const)), fmt(r$l1_grad),
              substr(r$l1_error %||% "", 1, 100)))
}
for (id in c("p116", "p053")) {
  run(id, "survey points: control", survey = TRUE)
  run(id, "survey points: z_1 row 2 reversed", 2L, survey = TRUE)
  run(id, "fixed points: control")
  run(id, "fixed points: z_1 row 2 reversed", 2L)
}
