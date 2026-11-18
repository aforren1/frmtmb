# Is frm() the slow step of the big_slope design? Timed alone.
source("dev/samplers-designs.R")
samplers_lib()
suppressPackageStartupMessages(library(frmtmb))
spec <- samplers_big("big_slope")
t0 <- proc.time()[[3]]
f <- frm(y ~ x + (1 + x | g), data = spec$data, family = bernoulli(),
         control = frmtmb_control(verbose = TRUE))
cat("frm seconds:", proc.time()[[3]] - t0, "\n")
print(f$opt$convergence); print(VarCorr(f))
