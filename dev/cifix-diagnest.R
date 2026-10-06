# The two fits that warn on Ubuntu and not on Windows at 0.68.0:
# test-aliased-grouping.R:138 and test-open-issues.R:52, both
# Reaction ~ Days + (1 | Subject/f) on sleepstudy. What does the SE
# check see on each BLAS?
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
data(sleepstudy, package = "lme4")
ss <- sleepstudy
ss$a <- factor(ss$Days %% 3)
for (f in c("a")) {
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(Reaction ~ Days + (1 | Subject/a)) + gaussian(), data = ss),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  cat("warnings:", length(w), substr(w, 1, 120), "\n")
  cat("par:", format(fit$opt$par, digits = 8), "\n")
  cat("evals:", fit$opt$evals, " conv:", fit$opt$convergence, "\n")
  an <- fit$cache$se_analysis
  cat("fit-time analysis tier:", if (is.null(an)) "none" else an$tier, "\n")
  h <- fit$cache$hessian_fixed
  if (!is.null(h)) {
    print(signif(h$H, 6))
    cat("noise E:\n"); print(signif(h$E, 3))
    u <- fit$par_units %||% rep(1, nrow(h$H))
    cat("row max:", signif(apply(abs(h$H * outer(u, u)), 1, max), 4), "\n")
    cat("10 * noise row max:",
        signif(10 * apply(h$E * outer(u, u), 1, max), 4), "\n")
  }
  print(frmtmb:::sdr_of(fit)$se_lost)
}
