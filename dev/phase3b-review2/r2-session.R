# Reviewer 2, item 3.3 (session =), my own construction, no fitting:
# (1) bandit2arm_dual: the two-session objective equals the sum of the
#     single-session objectives at random parameter vectors;
# (2) rows shuffled: the objective does not move;
# (3) a constant session column: the objective equals the one without
#     session =;
# (4) with (1 | id) the per-subject random effect spans both sessions:
#     the objective differs from the fit that splits subjects by
#     session (id:session), which gives each session its own effect.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.learn)})
mkd <- function(seed) frm_task_design("bandit2arm", n_subject = 5L, n_trial = 40L, seed = seed)
d1 <- mkd(11); d2 <- mkd(12); d1$session <- "a"; d2$session <- "b"
d <- rbind(d1, d2)
d$choice <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial, session = session), d,
                              pars = list(alpha = 0.3, tau = 2), seed = 3)[[1L]]$choice
f <- bf(choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1)
obj <- function(dd, fam) frm(f, family = fam, data = dd, dry_run = "objective")$obj
o2 <- obj(d, bandit2arm_dual(subject = id, trial = trial, session = session))
oa <- obj(d[d$session == "a", ], bandit2arm_dual(subject = id, trial = trial))
ob <- obj(d[d$session == "b", ], bandit2arm_dual(subject = id, trial = trial))
set.seed(4)
for (i in 1:3) {
  p <- o2$par + rnorm(length(o2$par), 0, 0.5)
  cat(sprintf("(1) dual: two-session %.12f  a + b %.12f  rel %.2e\n", o2$fn(p),
              oa$fn(p) + ob$fn(p), abs(o2$fn(p) - oa$fn(p) - ob$fn(p)) / abs(o2$fn(p))))
}
sh <- d[sample(nrow(d)), ]
os <- obj(sh, bandit2arm_dual(subject = id, trial = trial, session = session))
p <- o2$par + 0.2
cat(sprintf("(2) rows shuffled: %.12f vs %.12f  rel %.2e\n", os$fn(p), o2$fn(p),
            abs(os$fn(p) - o2$fn(p)) / abs(o2$fn(p))))
dc <- d1; dc$choice <- d$choice[d$session == "a"]
o_c <- obj(dc, bandit2arm_dual(subject = id, trial = trial, session = session))
o_n <- obj(dc, bandit2arm_dual(subject = id, trial = trial))
cat(sprintf("(3) constant session vs none: %.15g vs %.15g identical %s\n", o_c$fn(p), o_n$fn(p),
            identical(o_c$fn(p), o_n$fn(p))))
fr <- bf(choice | reward(pay1, pay2) ~ 1 + (1 | id), tau ~ 1)
d$idsess <- interaction(d$id, d$session)
r1 <- frm(fr, family = bandit2arm_delta(subject = id, trial = trial, session = session),
          data = d, dry_run = "objective")$obj
fr2 <- bf(choice | reward(pay1, pay2) ~ 1 + (1 | idsess), tau ~ 1)
r2 <- frm(fr2, family = bandit2arm_delta(subject = idsess, trial = trial),
          data = d, dry_run = "objective")$obj
cat(sprintf("(4) (1|id) over sessions: %.10f;  (1|id:session) as separate subjects: %.10f\n",
            r1$fn(r1$par), r2$fn(r2$par)))
cat("    random-effect lengths:", length(r1$env$random), "vs", length(r2$env$random), "\n")
