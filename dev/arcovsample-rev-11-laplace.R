# REVIEW script 11: posterior_predict() on a REAL laplace draws object
# came back all NaN in script 10. Is that this lane's change or does the
# reference build do the same?
#
#   Rscript dev/arcovsample-rev-11-laplace.R <lane|ref>

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]
LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
newx <- "arma_cond_resp" %in% getNamespaceExports("frmtmb")
cat("ARM ", arm, " newexports=", newx, "\n", sep = "")
stopifnot(identical(newx, identical(arm, "lane")))

set.seed(1212L)
dd <- data.frame(g = factor(rep(1:6, each = 5L)), t = rep(1:5, 6L))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(6L, 0, 0.5)[as.integer(dd$g)] +
  rnorm(nrow(dd), 0, 0.7)

one <- function(lab, form) {
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, family = gaussian(), data = dd, chains = 1,
               iter = 300, refresh = 0, seed = 3, laplace = TRUE)))
  set.seed(9)
  w <- NULL
  pp <- withCallingHandlers(
    tryCatch(posterior_predict(ds), error = identity),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning") })
  set.seed(9)
  ep <- tryCatch(posterior_epred(ds), error = function(e) e)
  cat("\n---- ", lab, "  laplace=", frmtmb.sample:::draws_is_laplace(ds),
      "  cols=", ncol(ds$draws), "\n", sep = "")
  cat("  draws columns: ", paste(colnames(ds$draws), collapse = " "),
      "\n", sep = "")
  if (inherits(pp, "condition")) {
    cat("  posterior_predict ERROR: ", conditionMessage(pp), "\n",
        sep = "")
  } else {
    cat("  posterior_predict: ", sum(!is.finite(pp)), " non-finite of ",
        length(pp), "\n", sep = "")
  }
  if (inherits(ep, "condition")) {
    cat("  posterior_epred ERROR: ", conditionMessage(ep), "\n", sep = "")
  } else {
    cat("  posterior_epred: ", sum(!is.finite(ep)), " non-finite of ",
        length(ep), "\n", sep = "")
  }
  cat("  distinct warnings (", length(unique(w)), "):\n", sep = "")
  for (s in head(unique(w), 4L)) cat("    ", s, "\n", sep = "")
}

one("ar(t, g) + (1 | g)", bf(y ~ x + ar(t, g) + (1 | g)))
one("(1 | g) only, no autocor", bf(y ~ x + (1 | g)))
one("ar(t, g) only, no random effect", bf(y ~ x + ar(t, g)))
cat("\nDONE\n")
