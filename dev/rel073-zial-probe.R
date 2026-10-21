# After the release review, finding 1: does zero_inflated_asym_laplace()
# with an observation-level effect have no maximum, as asym_laplace()
# has, so that the weights check may carry it in its table? Its log
# density on a row with y != 0 is log(1 - zi) plus asym_laplace's, linear
# on each side of the mode; a zero row adds log(zi) and no curvature.
#   Rscript dev/rel073-zial-probe.R <lib>   (log dev/rel073-log/zial-probe.txt)
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "\n")
run <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) {
    cat(sprintf("%-46s ERROR: %s\n", label,
                substr(gsub("[[:space:]]+", " ", r), 1, 110)))
  } else {
    cat(sprintf("%-46s logLik %10.4f  log sigma %8.3f  conv %d  warnings %d\n",
                label, as.numeric(logLik(r)), r$estimates$betad[1],
                r$opt$convergence, length(w)))
  }
}
for (seed in c(7, 8)) {
  set.seed(seed)
  n <- 80
  d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                  g = factor(rep(1:20, each = 4)))
  d$y <- ifelse(runif(n) < 0.2, 0, 1 + 0.5 * d$x + rnorm(n))
  d$w12 <- 1.2
  f <- zero_inflated_asym_laplace()
  cat("== seed", seed, "zeros", sum(d$y == 0), "of", n, "\n")
  run("(1 | obs), Laplace, unweighted",
      frm(y ~ x + (1 | obs), family = f, data = d))
  run("(1 | obs), quadrature, unweighted",
      frm(y ~ x + (1 | obs), family = f, data = d, quadrature = TRUE))
  run("(1 | obs), quadrature, weights 1.2",
      frm(y | weights(w12) ~ x + (1 | obs), family = f, data = d,
          quadrature = TRUE))
  run("(1 | g), 4 rows per level, Laplace",
      frm(y ~ x + (1 | g), family = f, data = d))
  run("no group term", frm(y ~ x, family = f, data = d))
}
cat("END\n")
