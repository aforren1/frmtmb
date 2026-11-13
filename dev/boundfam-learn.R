# Lane boundfam: frmtmb.learn test-engine.R's one-subject toy (25
# trials), whose tau the limit rule names at the end of its range. Is it?
#   Rscript dev/boundfam-learn.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.learn)
})
cat("lib", dirname(find.package("frmtmb")), "\n")
d <- frm_task_design("bandit2arm", n_subject = 1L, n_trial = 25L, seed = 44L)
d$choice <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
                              d, pars = list(alpha = 0.35, tau = 3),
                              seed = 44L)[[1L]]$choice
w <- character()
f <- withCallingHandlers(
  frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
      family = bandit2arm_delta(subject = id, trial = trial), data = d),
  warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
for (x in w) cat("W:", substr(x, 1, 300), "\n")
p <- f$opt$par
cat("code", f$opt$convergence, "par", signif(p, 6), "logLik",
    format(as.numeric(logLik(f)), digits = 12), "\n")
print(suppressWarnings(sqrt(diag(vcov(f, full = TRUE)))))
# the profile of tau, alpha refitted at each value
tr <- p[2] + c(-2, -1, 0, 2, 5, 10, 20)
for (t in tr) {
  o <- stats::optimize(function(al) f$obj$fn(c(al, t)), c(-30, 30))
  cat(sprintf("log tau %.3f: best -logLik %.10f at alpha %.4f\n", t,
              o$objective, o$minimum))
}
