# Review of lane nbrobust: where nbinom_cdf_log() turns NaN in the log
# shape, and whether a censored fit can reach it.
#   Rscript dev/reviews/nbrobust-review-cdfnan.R <base|lib>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ncl <- frmtmb:::nbinom_cdf_log
for (q in c(1, 3, 8, 20, 100)) for (lmu in c(-3, 0, 1.5, 4)) {
  ls <- seq(30, 60, by = 0.125)
  bad <- vapply(ls, function(s) {
    r <- ncl(q, lmu, s)
    !is.finite(r$lF) || !is.finite(r$lS)
  }, NA)
  first <- if (any(bad)) min(ls[bad]) else NA
  cat(sprintf("q %3d lmu %4g first non-finite log shape %s (of 30..60 by 1/8);",
              q, lmu, first),
      "finite again above it:", if (any(bad)) any(!bad[ls > first]) else NA,
      "\n")
}
cat("value at q 3 lmu 0, log shape 37..40:\n")
for (s in seq(37, 40, by = 0.5)) {
  r <- ncl(3, 0, s)
  cat(sprintf("  %5.2f lF %s lS %s\n", s, format(r$lF), format(r$lS)))
}
cat("stats::pnbinom(3, size = exp(40), mu = 1, log = TRUE):",
    stats::pnbinom(3, size = exp(40), mu = 1, log.p = TRUE), "\n")
# a censored fit, from log shape 38: does the objective read NaN
set.seed(2)
d <- data.frame(x = stats::rnorm(200))
d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
d$cc <- ifelse(d$y >= 5, "right", "none")
d$yc <- pmin(d$y, 5)
r <- tryCatch(withCallingHandlers(
  frm(yc | cens(cc) ~ x, family = negbinomial(), data = d,
      start = list(betad = 40)),
  warning = function(w) {
    cat("  WARNING:", substr(conditionMessage(w), 1, 160), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })
if (!is.null(r)) {
  cat(sprintf("cens from 40: code %d lsh %.4f logLik %.10f\n",
              r$opt$convergence, r$opt$par[["betad"]],
              as.numeric(logLik(r))))
}
