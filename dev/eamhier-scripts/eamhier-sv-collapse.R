# Lane eamhier: the replicate where `sv` collapsed, and what the
# post-fit surface says about it.
#
# Seed 20260935 of the single-level control (12,000 trials, one subject,
# sv = 0.4 in the simulator) returns sv = 4.4e-04 with a log-scale Wald
# interval of (-1163, +1147), convergence code 0, a positive definite
# Hessian and no NaN standard error. This script refits it and prints
# everything `diagnose()` has to say, because "nothing flags it" is a
# claim about the surface and has to be read off the surface.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-sv-collapse.R
source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

tr <- eamhier_truth
n <- 12000L
set.seed(20260935L)
cond <- rep(0:1, each = n / 2L)
d <- ddm_simulate(n, mu = tr$mu0 + tr$mu_cond * cond, bs = tr$bs,
                  ndt = tr$ndt, bias = 0.5, sv = tr$sv)
d$cond <- factor(cond, labels = c("a", "b"))

fit <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
           family = wiener(variability = "sv"), data = d, se = TRUE)
cat("logLik ", as.numeric(logLik(fit)), "\n")
print(unlist(fixef(fit)))
print(suppressWarnings(stats::confint(fit)))
cat("\n== diagnose(quiet = FALSE) ==\n")
dg <- frmtmb::diagnose(fit)
cat("\n== the fields, in full ==\n")
utils::str(dg)

# The same data with `sv` switched off, to say what the collapse cost.
f0 <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = wiener(), data = d, se = TRUE)
cat("\nlogLik with sv estimated ", as.numeric(logLik(fit)),
    "\nlogLik with no sv term   ", as.numeric(logLik(f0)),
    "\ndifference               ",
    as.numeric(logLik(fit)) - as.numeric(logLik(f0)), "\n")

# And at the truth, to say whether the optimizer stopped short or the
# likelihood really is flat there.
fx <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5,
             sv = 0.4),
          family = wiener(variability = "sv"), data = d, se = TRUE)
cat("logLik with sv held at 0.4", as.numeric(logLik(fx)), "\n")
