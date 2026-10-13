# Reviewer of lane optrap, final check: test-se-tier.R's healthy
# sleepstudy fit, now in units of 32 ms. Per scaling of Reaction: the
# unit, what the fit says, the SE check's lost set, and how far vcov()
# is from sdreport()'s cov.fixed (relative to the SEs).
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
for (k in c(1, 32, 64)) {
  ss <- sleepstudy
  ss$Reaction <- ss$Reaction / k
  n <- 0L
  f <- withCallingHandlers(frm(Reaction ~ Days + (Days | Subject), data = ss),
    warning = function(w) {n <<- n + 1L; invokeRestart("muffleWarning")},
    message = function(m) {n <<- n + 1L; invokeRestart("muffleMessage")})
  V <- vcov(f, full = TRUE)
  ref <- RTMB::sdreport(f$obj)$cov.fixed
  l4 <- lmer(Reaction ~ Days + (Days | Subject), ss, REML = FALSE)
  cat(sprintf("Reaction / %g | unit %s | said %d | lost %d | identical to sdreport %s | max rel %.3g | logLik - lmer %.3g\n",
              k, paste(signif(unique(f$par_units %||% 1), 3), collapse = ","),
              n, length(frmtmb:::sdr_of(f)$se_lost),
              identical(unname(V), unname(ref)),
              max(abs(V - ref) / sqrt(outer(diag(ref), diag(ref)))),
              as.numeric(logLik(f)) - as.numeric(logLik(l4))))
}
