# Lane bigx: frmtmb_control(autoscale = TRUE) on poisson
# y ~ (z + m) + (1 | g): the warnings in full, the code and the SEs.
#   BIGX_ARM=before|after Rscript dev/bigx-auto.R [m] [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
a <- commandArgs(TRUE)
m <- if (length(a)) as.numeric(a[1]) else 1e5
s <- if (length(a) > 1) as.integer(a[2]) else 1L
d <- bigx_data("poisson", m, s)
fa <- bigx_fit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
               control = frmtmb_control(autoscale = TRUE, verbose = 1))
cat("code", fa$f$opt$convergence, fa$f$opt$message, "\n")
cat("par", format(fa$f$opt$par, digits = 10), "\n")
cat("par_units", format(fa$f$par_units), "\n")
cat("SE", format(bigx_se(fa$f)), "\n")
w <- character(0)
withCallingHandlers(
  frm(bf(y ~ w + (1 | g)), family = poisson(), data = d,
      control = frmtmb_control(autoscale = TRUE)),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
cat("warnings:\n")
writeLines(strwrap(w, 78, prefix = "  "))
