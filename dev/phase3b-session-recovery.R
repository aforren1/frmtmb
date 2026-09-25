# Item 3.3: recovery and Wald coverage of session = at a realistic scale.
#
# Usage: Rscript dev/phase3b-session-recovery.R <seed_from> <seed_to>
# Design: 40 subjects, 2 sessions of 100 trials each (8000 rows),
# bandit2arm_delta. Each session is a fresh task, drawn from a fresh
# value store. alpha = plogis(qlogis(0.35) + u_id), u_id ~ N(0, 0.5),
# shared across a subject's two sessions; tau = 3.
# Two fits of the same data: WITH session = (the model that generated
# it) and WITHOUT (the store carried across the boundary), so the cost
# of ignoring the boundary is measured on the same draws.
# One RDS per seed in dev/phase3b-log/recov/session-<seed>.rds.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.learn)
})
args <- commandArgs(trailingOnly = TRUE)
seeds <- seq(as.integer(args[[1]]), as.integer(args[[2]]))
dir.create("dev/phase3b-log/recov", showWarnings = FALSE, recursive = TRUE)
NS <- 40L; NT <- 100L

for (seed in seeds) {
  t0 <- proc.time()[["elapsed"]]
  d1 <- frm_task_design("bandit2arm", n_subject = NS, n_trial = NT,
                        seed = seed)
  d2 <- frm_task_design("bandit2arm", n_subject = NS, n_trial = NT,
                        seed = seed + 100000L)
  d1$session <- 1L
  d2$session <- 2L
  d <- rbind(d1, d2)
  set.seed(seed)
  u <- rnorm(NS, 0, 0.5)
  al <- plogis(qlogis(0.35) + u[as.integer(factor(d$id))])
  # drawn with the subject:session pair as the learner, which is the
  # generative statement "each session starts from a fresh store"
  d$idsess <- interaction(d$id, d$session)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = idsess, trial = trial), d,
    pars = list(alpha = al, tau = 3), seed = seed)[[1]]$choice
  f <- bf(choice | reward(pay1, pay2) ~ 1 + (1 | id), tau ~ 1)
  one <- function(fam, dd) {
    tryCatch({
      fit <- frm(f, family = fam, data = dd)
      list(confint = confint(fit), fixef = unlist(fixef_by_dpar(fit)),
           varcorr = VarCorr(fit)$id$sd, loglik = as.numeric(logLik(fit)),
           conv = fit$opt$convergence, pdhess = fit$sdr$pdHess)
    }, error = function(e) list(error = conditionMessage(e)))
  }
  with_s <- one(bandit2arm_delta(subject = id, trial = trial,
                                 session = session), d)
  # without: the second session's trials numbered after the first's
  dn <- d
  dn$trial <- dn$trial + NT * (dn$session - 1L)
  without <- one(bandit2arm_delta(subject = id, trial = trial), dn)
  res <- list(seed = seed, with_session = with_s, without = without,
              elapsed = proc.time()[["elapsed"]] - t0,
              pkg = as.character(packageVersion("frmtmb.learn")))
  saveRDS(res, sprintf("dev/phase3b-log/recov/session-%d.rds", seed))
  cat("session", seed, round(res$elapsed, 1), "\n")
}
