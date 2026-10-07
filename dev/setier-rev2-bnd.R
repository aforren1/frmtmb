# Reviewer of lane setier, re-check: dev/setier-rev-bnd.R with the check_se = "stop" fit wrapped (it now stops at frm()).
# Reviewer of lane setier: the boundary message, once per fit, and what
# the post-fit surface says and returns on boundary fits.
#   Rscript dev/setier-rev-bnd.R lane|base
# A  sleepstudy (1 | Subject/a): conditions per call, fit-time check
# B  a deferred check: y ~ f (20 levels) + (1 | g), no g variance
# C  gaussian GAM, a smooth at its limit (test-se-tier.R's data): the
#    prediction surface
# D  cumulative() with (1 | g), no g variance: fitted() and predict()
# E  update() and frm(se = TRUE) on the A fit
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  merge = c("C:/Users/adf44/source/r/setier-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
grDevices::pdf(NULL)
conds <- function(lab, expr) {
  w <- character(); m <- character()
  v <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  bm <- grepl("^Boundary [(]singular[)] fit", m)
  cat(sprintf("  %-34s boundary msg %d | other msg %d | warnings %d %s%s\n",
              lab, sum(bm), sum(!bm), length(w),
              if (length(w)) paste0("[", substr(w[1], 1, 70), "]") else "",
              if (inherits(v, "err")) paste0(" ERROR: ", substr(v, 1, 80))
              else ""))
  invisible(v)
}
nanc <- function(x) sum(is.nan(as.numeric(x)))

cat("\n== A: sleepstudy (1 | Subject/a)\n")
d <- lme4::sleepstudy
d$a <- factor(d$Days %% 3)
fA <- conds("frm()", frm(Reaction ~ Days + (1 | Subject/a), data = d))
cat("  deferred:", !is.null(fA$cache$se_deferred), "\n")
conds("summary()", summary(fA))
conds("print(summary())", capture.output(print(summary(fA))))
conds("vcov()", vcov(fA))
conds("fixef()", fixef(fA))
conds("VarCorr()", VarCorr(fA))
conds("confint()", confint(fA))
conds("confint_varcorr()", confint_varcorr(fA))
r <- conds("ranef(condVar = TRUE)", ranef(fA, condVar = TRUE))
p <- conds("predict(ndraws = 100)", predict(fA, ndraws = 100))
cat("  predict NaN Est.Error:", nanc(p[, "Est.Error"]), "of", nrow(p), "\n")
p <- conds("fitted()", fitted(fA))
cat("  fitted NaN Est.Error:", nanc(p[, "Est.Error"]), "\n")
p <- conds("frm_linpred(se.fit)", frm_linpred(fA, se.fit = TRUE))
cat("  linpred NaN se:", nanc(p$se.fit), "\n")
ce <- conds("conditional_effects()", conditional_effects(fA))
cat("  ce NaN se__:", nanc(ce[[1]]$se__), "\n")
conds("diagnose()", diagnose(fA, quiet = TRUE))
conds("hypothesis(Days > 0)", hypothesis(fA, "Days > 0"))
conds("emmeans", if (requireNamespace("emmeans", quietly = TRUE))
  emmeans::emmeans(fA, ~ 1) else NULL)
cat("\n== E: update(), se = TRUE, a second frm()\n")
conds("update(fA)", update(fA))
conds("frm(se = TRUE)", frm(Reaction ~ Days + (1 | Subject/a), data = d,
                            se = TRUE))

cat("\n== B: deferred check, y ~ f + (1 | g), 20-level f, no g variance\n")
set.seed(3)
dB <- data.frame(f = factor(rep(1:20, 6)), g = factor(rep(1:8, 15)))
dB$y <- rnorm(20)[dB$f] + rnorm(120)
fB <- conds("frm()", frm(y ~ f + (1 | g), data = dB))
cat("  outer pars", length(fB$opt$par), "evals", fB$opt$evals, "deferred:",
    !is.null(fB$cache$se_deferred), "\n")
conds("fixef() (first SE use)", fixef(fB))
conds("summary()", summary(fB))
conds("VarCorr()", VarCorr(fB))
conds("vcov()", vcov(fB))
fB2 <- suppressMessages(frm(y ~ f + (1 | g), data = dB))
conds("suppressMessages(frm()) then VarCorr()", VarCorr(fB2))
conds("check_se = stop: frm()", frm(y ~ f + (1 | g), data = dB,
           control = frmtmb_control(check_se = "stop")))

cat("\n== C: gaussian GAM, s(x) + s(z)\n")
set.seed(3)
dC <- data.frame(x = runif(200), z = runif(200))
dC$y <- 2 * dC$x + sin(2 * pi * dC$z) + rnorm(200, 0, 0.3)
fC <- conds("frm()", frm(y ~ s(x) + s(z), data = dC))
cat("  lost:", paste(names(ns$sdr_of(fC)$se_lost), ns$sdr_of(fC)$se_lost),
    "\n")
p <- conds("frm_linpred(se.fit)", frm_linpred(fC, se.fit = TRUE))
cat("  linpred NaN se:", nanc(p$se.fit), "of", length(p$se.fit), "\n")
p <- conds("fitted()", fitted(fC))
cat("  fitted NaN Est.Error:", nanc(p[, "Est.Error"]), "\n")
p <- conds("predict(ndraws = 100)", predict(fC, ndraws = 100))
ce <- conds("conditional_effects()", conditional_effects(fC))
cat("  ce NaN se__:", paste(vapply(ce, function(x) nanc(x$se__), 0),
                            collapse = " "), "\n")
s <- conds("summary()", summary(fC))
conds("smooth_terms / VarCorr()", VarCorr(fC))
if (requireNamespace("mgcv", quietly = TRUE)) {
  g <- mgcv::gam(y ~ s(x) + s(z), data = dC, method = "ML")
  cat("  mgcv sp:", signif(g$sp, 3), " edf:", signif(summary(g)$edf, 3),
      "\n")
  pg <- predict(g, se.fit = TRUE)
  pf <- frm_linpred(fC, se.fit = TRUE)
  cat("  linpred SE frmtmb / mgcv (ML), median and range:",
      signif(median(pf$se.fit / pg$se.fit), 3),
      signif(range(pf$se.fit / pg$se.fit), 3), "\n")
}

cat("\n== D: cumulative() with (1 | g), no g variance\n")
set.seed(9)
dD <- data.frame(x = rnorm(150), g = factor(rep(1:10, 15)))
dD$y <- cut(0.8 * dD$x + rlogis(150), c(-Inf, -1, 0.5, 1.5, Inf),
            labels = FALSE)
fD <- conds("frm()", frm(y ~ x + (1 | g), family = cumulative(),
                         data = dD))
cat("  lost:", paste(names(ns$sdr_of(fD)$se_lost), ns$sdr_of(fD)$se_lost),
    "\n")
p <- conds("fitted()", fitted(fD))
cat("  fitted NaN Est.Error:", nanc(p[, "Est.Error", ]), "of",
    length(p[, "Est.Error", ]), "\n")
p <- conds("fitted(re_formula = NA)", fitted(fD, re_formula = NA))
cat("  fitted NaN Est.Error:", nanc(p[, "Est.Error", ]), "\n")
p <- conds("predict(ndraws = 100)", predict(fD, ndraws = 100))
ce <- conds("conditional_effects()", conditional_effects(fD))
