# Lane optstart: where the re1000 false convergence comes from
# (seeds 6 and 14 of dev/optrap-rev2-b1.R re1000): the stages of the fit
# with verbose = TRUE, the first run's code and gradient in natural
# units, and what the restart does.
#   Rscript dev/optstart-re1000-trace.R <lib or ""> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else c(6, 14)
cat("lib", dirname(find.package("frmtmb")), "\n")
for (s in seeds) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 1000)[d$g] + rnorm(100)
  cat("=== seed", s, "\n")
  f0 <- suppressWarnings(frm(y ~ x + (1 | g), data = d,
                             control = frmtmb_control(restarts = 0)))
  u <- f0$par_units %||% rep(1, length(f0$opt$par))
  cat("restarts = 0: code", f0$opt$convergence, f0$opt$message,
      "| iter", f0$opt$iterations, "| max|grad * unit|",
      max(abs(f0$obj$gr(f0$opt$par) * u)), "| units",
      paste(signif(u, 3), collapse = ","), "\n")
  cat("   start (par_template):", format(unlist(par_template(
    y ~ x + (1 | g), data = d)[c("beta", "betad", "theta")]), digits = 5),
    "\n")
  f1 <- withCallingHandlers(frm(y ~ x + (1 | g), data = d,
                                control = frmtmb_control(verbose = TRUE)),
    warning = function(w) {
      cat("   WARNING:", substr(conditionMessage(w), 1, 90), "\n")
      invokeRestart("muffleWarning")
    })
  cat("default: code", f1$opt$convergence, "| logLik diff",
      as.numeric(logLik(f1)) - as.numeric(logLik(f0)), "\n")
}
