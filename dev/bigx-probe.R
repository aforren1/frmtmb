# Lane bigx, step 1: what goes wrong on poisson y ~ (z + m) + (1 | g).
# For each offset m and data seed: the fit on w = z + m against the same
# model on z (an exact reparametrization), the stop's intermediates
# (code, message, iterations, gradient, Newton polish, design kappa,
# autoscale), and autoscale = TRUE.
#   BIGX_ARM=before|after Rscript dev/bigx-probe.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
ns <- asNamespace("frmtmb")
for (m in c(1e3, 1e4, 1e5, 1e6, 1e7)) for (s in 1:3) {
  d <- bigx_data("poisson", m, s)
  fc <- bigx_fit(bf(y ~ z + (1 | g)), family = poisson(), data = d)
  fw <- bigx_fit(bf(y ~ w + (1 | g)), family = poisson(), data = d)
  fa <- bigx_fit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
                 control = frmtmb_control(autoscale = TRUE))
  sec <- bigx_se(fc$f)
  line <- function(lbl, x) {
    if (x$err) {
      cat(sprintf("  %-4s STOP %s\n", lbl, substr(conditionMessage(x$f), 1,
                                                  100)))
      return(invisible())
    }
    f <- x$f
    b <- f$opt$par
    g <- tryCatch(as.numeric(f$obj$gr(b)), error = function(e) NA)
    se <- bigx_se(f)
    cat(sprintf(paste0("  %-4s slope z %+9.3f theta z %+7.3f dll %+.3e ",
                       "SE ratio %.4f code %d it %s |g| %.2e polish %s ",
                       "msg '%s' warn %d [%s]\n"),
                lbl, (b[[2]] - fc$f$opt$par[[2]]) / sec[[2]],
                (b[[3]] - fc$f$opt$par[[3]]) / sec[[3]],
                as.numeric(logLik(f)) - as.numeric(logLik(fc$f)),
                se[[2]] / sec[[2]], f$opt$convergence,
                format(f$opt$iterations), max(abs(g)),
                if (is.null(f$opt$polish)) "-" else
                  sprintf("%d/%.2e", f$opt$polish[["steps"]],
                          f$opt$polish[["gain"]]),
                substr(f$opt$message, 1, 30), length(x$w),
                paste(unique(x$wc), collapse = ",")))
  }
  fr <- if (!fw$err) fw$f$frame
  kap <- if (!is.null(fr)) ns$predictor_design_kappa(fr) else NA
  cat(sprintf("m %g seed %d kappa %.3e | centered slope %.5f SE %.5f\n",
              m, s, kap, fc$f$opt$par[[2]], sec[[2]]))
  line("raw", fw)
  line("auto", fa)
  if (length(fw$w)) cat("     raw warns:", paste(fw$w, collapse = " || "), "\n")
}
