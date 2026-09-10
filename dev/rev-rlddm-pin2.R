# REVIEW round 2, attack 1: the TWO justifications the refusal now
# rests on, since they are the whole case for a breaking change.
#
#   Rscript dev/rev-rlddm-pin2.R <lib>
#
# Seed 4242, the design dev/rlddm-scripts/rlddm-pin.R uses.
#
# A. `bf(ndt = 0.3)` with `rlddm(max_ndt = 0.26)`. Run through a REAL
#    frm() rather than dry_run, because the claim is about what a user
#    SEES, and a dry run does not optimize.
# B. Consistency: `wiener()` refuses `bf(ndt = 0.15)` and
#    `wiener(max_ndt = 0.2)` accepts it, on BOTH builds. If that holds,
#    leaving rlddm() different would make one family of five behave
#    unlike the other four.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
cat("lib  :", lib, "\n")
cat("learn:", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n")
cat("eam  :", format(packageVersion("frmtmb.eam")), "at",
    dirname(system.file(package = "frmtmb.eam")), "\n\n")

run <- function(lab, e) {
  out <- withCallingHandlers(
    tryCatch({
      v <- e
      paste("RETURNED, logLik =",
            format(as.numeric(stats::logLik(v)), digits = 9))
    }, error = function(err)
      paste("ERROR:", gsub("[\r\n]+", " ", conditionMessage(err)))),
    warning = function(w) invokeRestart("muffleWarning"))
  cat(sprintf("%-42s %s\n", lab, substr(out, 1, 190)))
  invisible(out)
}

set.seed(4242)
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 60,
                     seed = 4242)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5),
                       seed = 4242)[[1L]]
cat("rlddm design, min(rt):", format(min(s$rt), digits = 9), "\n")

f3 <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1, bs ~ 1,
         ndt = 0.3, bias = 0.5)
run("A. bf(ndt = 0.3), rlddm(max_ndt = 0.26), frm()",
    frm(f3, family = rlddm(subject = id, trial = trial,
                           max_ndt = 0.26), data = s))
f2 <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1, bs ~ 1,
         ndt = 0.2, bias = 0.5)
run("   the same at bf(ndt = 0.2), which is in range",
    frm(f2, family = rlddm(subject = id, trial = trial,
                           max_ndt = 0.26), data = s))
run("   bf(ndt = 0.2), no max_ndt",
    frm(f2, family = rlddm(subject = id, trial = trial), data = s))

cat("\nB. the four eam families, same question\n")
set.seed(11)
w <- ddm_simulate(400, mu = 1.2, bs = 1.5, ndt = 0.2)
cat("   wiener design, min(rt):", format(min(w$rt), digits = 9), "\n")
wf <- bf(rt | dec(upper) ~ 1, bs ~ 1, ndt = 0.15, bias = 0.5)
run("   bf(ndt = 0.15), wiener()",
    frm(wf, family = wiener(), data = w))
run("   bf(ndt = 0.15), wiener(max_ndt = 0.2)",
    frm(wf, family = wiener(max_ndt = 0.2), data = w))
wf3 <- bf(rt | dec(upper) ~ 1, bs ~ 1, ndt = 0.25, bias = 0.5)
run("   bf(ndt = 0.25), wiener(max_ndt = 0.2), out of range",
    frm(wf3, family = wiener(max_ndt = 0.2), data = w))
run("   bf(ndt = 0.15), rdm(2)",
    frm(bf(rt | dec(upper) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt = 0.15),
        family = rdm(2), data = w))
run("   bf(ndt = 0.15), rdm(2, max_ndt = 0.2)",
    frm(bf(rt | dec(upper) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt = 0.15),
        family = rdm(2, max_ndt = 0.2), data = w))
