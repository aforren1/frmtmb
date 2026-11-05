# The learn fit's profile returned NA for the lower bound with no warning. Does
# the Laplace fit do the same, and is the likelihood at sd -> 0 above the cutoff?
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.learn)
})
ln_small <- function(ns, nt, seed, sd_u) {
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt, seed = seed)
  set.seed(seed)
  al <- stats::plogis(stats::qlogis(0.4) +
                        stats::rnorm(ns, 0, sd_u)[as.integer(d$id)])
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = al, tau = 3), seed = seed)[[1L]]$choice
  d
}
d <- ln_small(20L, 60L, 63L, 0.8)
f <- bf(choice | reward(pay1, pay2) ~ 1 + (1 | id), tau ~ 1)
fam <- bandit2arm_delta(subject = id, trial = trial)
lap <- suppressWarnings(frm(f, family = fam, data = d))
imp <- suppressWarnings(frm(f, family = fam, data = d, importance = 400L))
cat("Laplace profile:\n")
print(suppressWarnings(confint(lap, parm = "theta_1", method = "profile")))
cat("importance profile:\n")
print(suppressWarnings(confint(imp, parm = "theta_1", method = "profile")))
# the corrected objective against theta, reoptimized over the others on the tape
cut <- imp$opt$objective + 0.5 * qchisq(0.95, 1)
cat("objective at the fit", imp$opt$objective, "cutoff", cut, "\n")
for (th in c(-1, -1.5, -2, -2.5, -3, -4, -6)) {
  o <- nlminb(imp$opt$par[1:2], function(b) imp$obj$fn(c(b, th)),
              function(b) imp$obj$gr(c(b, th))[1:2])
  l <- nlminb(lap$opt$par[1:2], function(b) lap$obj$fn(c(b, th)),
              function(b) lap$obj$gr(c(b, th))[1:2])
  cat(sprintf("theta %5.1f importance tape %.4f (above cutoff %s)  Laplace %.4f\n", th,
              o$objective, o$objective > cut, l$objective))
}
