# Reviewer scratch (fix-round check of lane translator): a half-collapsed
# block. p087 (gr_cov x mmc()) has z_1 = matrix[2, levels] with SDs
# 0.463 / 1.12e-06 (dev/translator-audit-out/p087.tsv). The audit reverses
# the WHOLE block, which row 1 carries, and calls it seen. Here only row 2
# (the mmc() slope) is reversed, after the translator, at the survey's
# points and at the lane's fixed points (scale 0.3). Records nothing.
#
#   Rscript dev/review-translator-p087.R
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
id <- "p087"
row <- smp[smp$id == id, ]
key <- paste(row$feature_a, row$feature_b, sep = " x ")
run <- function(label, rev_row = NULL, points = FALSE) {
  e <- new.env(parent = asNamespace("frmtmb"))
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
    sys.source(f, envir = e)
  }
  if (points) {
    mk <- e$eq_specs[[key]]
    e$eq_specs[[key]] <- function() {
      s <- mk()
      s$l1_points <- list(scale = 0.3, why = "review")
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
  cat(sprintf("%-40s verdict %-12s L1 ok %-5s const %s grad %s err %s\n",
              label, r$verdict %||% "", format(r$l1_ok %||% NA),
              fmt(c(r$l1_const, r$l1_extra_const)), fmt(r$l1_grad),
              substr(r$l1_error %||% "", 1, 100)))
}
run("survey points: control")
run("survey points: z_1 row 2 reversed", 2L)
run("survey points: z_1 row 1 reversed", 1L)
run("fixed points: control", points = TRUE)
run("fixed points: z_1 row 2 reversed", 2L, points = TRUE)
