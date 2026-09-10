# REVIEW: the refusal-first ordering in rlddm()'s family_finalize, with
# the guarded thing ABSENT: a family whose bound is already SETTLED at
# construction would skip the check entirely if the carried bound were
# read first.
.libPaths(c("C:/Users/adf44/source/r/rev-rlddm-lib","C:/Users/adf44/source/r/pinlib","C:/Users/adf44/source/r/rellib-0552","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb);library(frmtmb.eam);library(frmtmb.learn)})
set.seed(4242)
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 60, seed = 4242)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
  pars = list(alpha = 0.4, drift = 3, bs = 1.6, ndt = 0.2, bias = 0.5),
  seed = 4242)[[1L]]
cat("min(rt):", format(min(s$rt), digits = 9), "\n")
f <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)
for (mx in c(5, 0.26)) {
  r <- tryCatch({ frm(f, family = rlddm(subject = id, trial = trial, max_ndt = mx), data = s, dry_run = "objective"); "ACCEPTED" },
                error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 90)))
  cat(sprintf("max_ndt = %-6s %s\n", mx, r))
}
# and the second finalize, run by hand on fewer rows, on a settled bound
o <- frm(f, family = rlddm(subject = id, trial = trial, max_ndt = 0.26), data = s, dry_run = "objective")
fam <- frmtmb::single_response(o)[["family"]]
drop <- s[-which.min(s$rt), ]
again <- fam[["family_finalize"]](fam, drop$rt, list())
cat("bound kept across a refit on fewer rows:", identical(again[["ndt_bound"]][["ub"]], fam[["ndt_bound"]][["ub"]]),
    " ub =", format(again[["ndt_bound"]][["ub"]], digits = 9), "\n")