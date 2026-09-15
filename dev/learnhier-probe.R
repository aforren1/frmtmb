# Lane `learnhier`: what the correlated block looks like from outside.
#
# Answers, on a toy design so it costs seconds:
#   1  the ORDER of the block's rows, which is the order the Stan
#      program must declare its `u` matrix in;
#   2  whether the `b` entries of the parameter vector are level-major
#      in that same order;
#   3  what `bf(bias = 0.5)` leaves in the parameter vector;
#   4  which rows confint() gives for the standard deviations and the
#      correlations, which is what a coverage table reads.
source("dev/learnhier-env.R")
suppressPackageStartupMessages(library(frmtmb.learn))
lane_env_report()

set.seed(11)
ns <- 8L
d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = 60, seed = 11)
ua <- stats::rnorm(ns, 0, 0.6)
ut <- stats::rnorm(ns, 0, 0.3)
i <- as.integer(d$id)
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = id, trial = trial), d,
  pars = list(alpha = stats::plogis(stats::qlogis(0.35) + ua[i]),
              tau = exp(log(3) + ut[i])), seed = 11)[[1L]]$choice

fit <- frmtmb::frm(
  frmtmb::bf(choice | reward(pay1, pay2) ~ 1 + (1 | p | id),
             tau ~ 1 + (1 | p | id)),
  family = bandit2arm_delta(subject = id, trial = trial), data = d)

cat("\n== bandit2arm_delta, (1 | p | id) on both ==\n")
V <- frmtmb::VarCorr(fit)
cat("VarCorr names: ", paste(names(V), collapse = ", "), "\n")
print(V[[1L]])
cat("rownames: ", paste(rownames(V[[1L]]), collapse = " | "), "\n")

p <- fit$obj$env$last.par.best
cat("par name table: \n")
print(table(names(p)))
cat("first 8 b entries:\n")
print(round(unname(p[names(p) == "b"])[1:8], 4))
cat("re_blocks:\n")
str(fit$frame$re_blocks, max.level = 2L)

cat("\nconfint rows:\n")
print(rownames(suppressWarnings(stats::confint(fit))))
cat("\nranef structure:\n")
r <- frmtmb::ranef(fit)
str(r, max.level = 2L)
print(utils::head(r[[1L]], 3L))
