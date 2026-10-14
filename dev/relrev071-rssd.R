# Release reviewer of 0.71.0: the false-alarm battery's rs_sd0.02 seed
# 4 (and any other seed given), the one new convergence warning of the
# merged build. Same data stream as dev/illse-sf-rev3-falarm.R's
# "mixed" group (the ri_sd fit runs first, as there).
#   Rscript dev/relrev071-rssd.R <lib> <seeds> [sdg]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "| R", R.home(), "\n")
seeds <- eval(parse(text = a[2]))
sdw <- if (length(a) > 2) as.numeric(a[3]) else 0.02
cap <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
  })
  list(v = v, w = w, m = m)
}
for (s in seeds) {
  set.seed(s)
  for (sdg in c(0.02, 0.05, 0.1)) {
    d <- data.frame(x = rnorm(300), g = factor(rep(1:30, 10)))
    d$y <- 1 + 0.5 * d$x + rnorm(30, 0, sdg)[d$g] + rnorm(300)
    invisible(cap(frm(y ~ x + (1 | g), data = d)))
    d$y2 <- 1 + 0.5 * d$x + rnorm(30, 0, 0.5)[d$g] +
      rnorm(30, 0, sdg)[d$g] * d$x + rnorm(300)
    if (sdg == sdw) {
      r <- cap(frm(y2 ~ x + (1 + x | g), data = d))
      f <- r$v
      m <- suppressMessages(lme4::lmer(y2 ~ x + (1 + x | g), data = d,
                                       REML = FALSE))
      cat(sprintf("\nseed %d sd %g: code %d (%s) iters %s evals %s\n", s,
                  sdg, f$opt$convergence, f$opt$message,
                  f$opt$iterations, paste(f$opt$evaluations, collapse = ",")))
      cat(sprintf("  logLik frm %.10f  lmer %.10f  diff %.3g\n",
                  as.numeric(logLik(f)), as.numeric(logLik(m)),
                  as.numeric(logLik(f)) - as.numeric(logLik(m))))
      cat("  theta", format(f$estimates$theta, digits = 8), "\n")
      cat("  lmer sd", format(as.data.frame(lme4::VarCorr(m))$sdcor,
                              digits = 6), "\n")
      if (!is.null(f$opt$restart)) {
        cat("  restart:", paste(names(f$opt$restart), collapse = ","), "\n")
        str(f$opt$restart[setdiff(names(f$opt$restart), "par")])
      }
      g <- f$obj$gr(f$opt$par)
      cat("  grad", format(drop(g), digits = 3), "\n")
      cat("  warnings:", if (length(r$w)) paste0("\n    ", r$w) else "none",
          "\n")
      cat("  messages:", if (length(r$m)) paste0("\n    ", r$m) else "none",
          "\n")
      se <- sqrt(diag(as.matrix(suppressWarnings(vcov(f, full = TRUE)))))
      cat("  SE", format(se, digits = 5), "\n")
      cat("  lmer SE", format(sqrt(diag(as.matrix(vcov(m)))), digits = 5),
          "\n")
      l <- ns$sdr_of(f)$se_lost
      cat("  lost:", if (length(l)) paste(names(l), l) else "-", "\n")
    }
  }
}
