# zero_inflated_negbinomial on test-nl-family.R's data with the zero
# draws of nl_fam_data() and the counts from rnbinom(size = s): which
# distributional parameter sits at the boundary in the intercept-only
# fit, per size.
#
#   Rscript dev/ci-zinbprobe.R <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- parse("tests/testthat/test-nl-family.R")
env <- new.env(parent = asNamespace("frmtmb"))
eval(src[[1L]], env)  # nl_fam_data
for (s in c(2, 3, 5, 8)) {
  z <- env$nl_fam_data()
  y <- ifelse(z$zero, 0L, stats::rnbinom(z$n, size = s, mu = exp(z$eta)))
  d <- data.frame(x = z$x, y = y)
  for (rhs in c("1", "x")) {
    w <- character()
    fit <- withCallingHandlers(
      suppressMessages(frm(bf(stats::as.formula(paste("y ~", rhs))),
                           family = zero_inflated_negbinomial(), data = d)),
      warning = function(cnd) {
        w <<- c(w, substr(conditionMessage(cnd), 1, 50))
        invokeRestart("muffleWarning")
      })
    e <- fit$estimates
    cat(sprintf("size %g, y ~ %s: shape %s zi %s ll %.8f %s\n", s, rhs,
                paste(signif(unlist(e[["betad"]])[1], 4)),
                paste(signif(unlist(e[["betad"]])[-1], 4), collapse = ","),
                as.numeric(logLik(fit)), paste(w, collapse = " | ")))
  }
}
