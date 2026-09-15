# Lane `learnhier`: does the toolchain RUN, rather than merely resolve?
#
#   Rscript dev/learnhier-smoke.R
#
# A version string is not evidence that a library is intact. This
# machine's R library has been destroyed five times, and the failure
# looks like a hollow package directory: `packageVersion()` still
# answers from the DESCRIPTION while the namespace will not load and a
# fit dies. So this fits both families the lane depends on, at toy
# sizes, and prints their log-likelihoods to ten digits so a rerun can
# be compared rather than eyeballed.
source("dev/learnhier-env.R")
lane_env_report()
suppressPackageStartupMessages(library(frmtmb.learn))

d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 40, seed = 1)
s <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, tau = 3), seed = 1)[[1L]]
f1 <- frmtmb::frm(
  frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
  family = bandit2arm_delta(subject = id, trial = trial), data = s)
cat("SMOKE bandit2arm_delta logLik ",
    format(as.numeric(stats::logLik(f1)), digits = 10),
    " conv ", f1$opt$convergence, "\n", sep = "")

r <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5),
                       seed = 1)[[1L]]
f2 <- frmtmb::frm(
  frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
             drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
  family = rlddm(subject = id, trial = trial), data = r)
cat("SMOKE rlddm logLik ",
    format(as.numeric(stats::logLik(f2)), digits = 10),
    " conv ", f2$opt$convergence, "\n", sep = "")
