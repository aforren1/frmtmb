# Lane setier: a boundary fit's downstream calls say nothing more than
# the fit's boundary message (sleepstudy, (1 | Subject/a)).
#   Rscript dev/setier-predcheck.R <lib or "base">
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
d <- lme4::sleepstudy
d$a <- factor(d$Days %% 3)
fit <- suppressMessages(suppressWarnings(
  frm(Reaction ~ Days + (1 | Subject/a), data = d)))
for (nm in c("predict", "VarCorr", "confint_varcorr", "confint",
             "conditional_effects", "fitted", "summary")) {
  w <- character()
  v <- withCallingHandlers(switch(nm,
    predict = predict(fit, ndraws = 200),
    VarCorr = VarCorr(fit), confint_varcorr = confint_varcorr(fit),
    confint = confint(fit),
    conditional_effects = conditional_effects(fit, "Days"),
    fitted = fitted(fit), summary = summary(fit)),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  cat(nm, "warnings", length(w), if (length(w)) substr(w[1], 1, 120), "\n")
}
