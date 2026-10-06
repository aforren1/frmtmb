# Reviewer: a REML fit (joint precision route) with a variance at 0.
# Does the SE warning's "the other standard errors are reported" hold,
# and how many warnings does the user see across frm() + summary()?
#   Rscript dev/nanse-rev-reml.R base|lane|merge
args <- commandArgs(trailingOnly = TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
n_se <- 0L
for (s in 1:30) {
  set.seed(200 + s)
  dd <- data.frame(g = factor(rep(1:8, each = 6)), x = rnorm(48))
  dd$y <- 1 + 0.5 * dd$x + rnorm(48)
  for (reml in c(FALSE, TRUE)) {
    r <- cap(frm(bf(y ~ x + (1 | g)), data = dd, REML = reml))
    v <- cap(vcov(r$v))
    se <- sqrt(diag(v$v))
    any_se <- any(grepl("^Standard errors are not available", r$w))
    if (any_se || length(v$w) || any(!is.finite(se))) {
      cat(sprintf("seed %d REML %s: code %d sd_g %.2g fixef SE %s\n",
                  200 + s, reml, r$v$opt$convergence,
                  exp(r$v$opt$par[["theta"]] %||% NA),
                  paste(signif(se, 3), collapse = " ")))
      for (x in r$w) cat("   frm():", substr(x, 1, 150), "\n")
      for (x in v$w) cat("   vcov():", substr(x, 1, 150), "\n")
    }
  }
}
