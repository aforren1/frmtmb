# Fix round 1: do exgaussian and asym_laplace with an observation-level
# intercept run away, and under which weights and which integration?
#   Rscript dev/fuzzfix-fr1-families.R <lib|base>
# Run under dev/fuzzfix-rev-capped.ps1 (small: n = 80, one block).
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261020)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                g = factor(rep(1:20, each = 4)))
eta <- 0.3 + 0.4 * d$x + rnorm(n, 0, 0.3)
d$ye <- eta + rnorm(n, 0, 0.4) + rexp(n, 1 / 0.5)
d$ya <- eta + rnorm(n, 0, 0.5)
run <- function(lbl, fo, fam, dd, quad = FALSE) {
  w <- character()
  f <- withCallingHandlers(
    tryCatch(frm(fo, family = fam, data = dd, quadrature = quad),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (inherits(f, "error")) {
    cat(sprintf("ARM %-34s ERROR %s\n", lbl,
                substr(gsub("\n", " ", conditionMessage(f)), 1, 100)))
  } else {
    cat(sprintf("ARM %-34s logLik %.6g betad %s conv %d | warn %d %s\n",
                lbl, as.numeric(logLik(f)),
                paste(format(f$estimates$betad, digits = 4), collapse = ","),
                f$opt$convergence, length(w),
                substr(paste(unique(w), collapse = " ; "), 1, 60)))
  }
}
for (wv in c(NA, 0.9, 1.2)) {
  dd <- d; dd$w <- if (is.na(wv)) 1 else wv
  lhs <- if (is.na(wv)) "" else " | weights(w)"
  tag <- if (is.na(wv)) "unweighted" else paste0("w=", wv)
  run(paste("exgaussian obs", tag),
      stats::as.formula(paste0("ye", lhs, " ~ x + (1 | obs)")),
      exgaussian(), dd)
  run(paste("exgaussian group4", tag),
      stats::as.formula(paste0("ye", lhs, " ~ x + (1 | g)")),
      exgaussian(), dd)
  run(paste("asym_laplace obs", tag),
      stats::as.formula(paste0("ya", lhs, " ~ x + (1 | obs)")),
      asym_laplace(), dd)
  run(paste("asym_laplace obs quad", tag),
      stats::as.formula(paste0("ya", lhs, " ~ x + (1 | obs)")),
      asym_laplace(), dd, quad = TRUE)
  run(paste("asym_laplace group4", tag),
      stats::as.formula(paste0("ya", lhs, " ~ x + (1 | g)")),
      asym_laplace(), dd)
}
cat("END\n")
