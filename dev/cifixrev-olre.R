# Reviewer: test-se-check.R's OLRE fixture, where sigma AND the id sd
# lose their standard errors on a fit with random effects, under ML and
# REML. Do predictions of sigma get NaN and the warning (lane nanse B2),
# and does REML's vcov() agree with the prediction covariance?
# Usage: Rscript dev/cifixrev-olre.R <lib|base>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
f6 <- function(v) format(signif(v, 6))
warns <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}
set.seed(5)
d2 <- data.frame(id = factor(1:80), x = rnorm(80))
d2$y <- 1 + 0.5 * d2$x + rnorm(80)
nd <- d2[1:3, ]
for (reml in c(FALSE, TRUE)) {
  fw <- warns(frm(bf(y ~ x + (1 | id)), family = gaussian(), data = d2,
                  REML = reml,
                  control = frmtmb_control(check_olre = "ignore")))
  fit <- fw$value
  cat("\n== REML", reml, "lost:",
      paste(names(ns$sdr_of(fit)$se_lost), ns$sdr_of(fit)$se_lost), "\n")
  for (w in fw$warnings) cat("  W:", substr(w, 1, 400), "\n")
  for (dp in c("mu", "sigma")) {
    r <- warns(frm_linpred(fit, newdata = nd, dpar = dp, se.fit = TRUE,
                           re_formula = NA))
    cat("frm_linpred dpar", dp, "se:", f6(r$value$se.fit), "warnings:",
        length(r$warnings), "\n")
    for (w in r$warnings) cat("  W:", substr(w, 1, 160), "\n")
  }
  r <- warns(fitted(fit, newdata = nd, dpar = "sigma", re_formula = NA))
  cat("fitted dpar sigma Est.Error:", f6(r$value[, "Est.Error"]),
      "warnings:", length(r$warnings), "\n")
  v <- warns(vcov(fit))
  cat("vcov diag:", f6(diag(v$value)), "warnings:", length(v$warnings), "\n")
  jc <- ns$get_joint_cov(fit)
  bpos <- which(jc$names == "beta")
  cat("get_joint_cov beta diag:", f6(diag(jc$V)[bpos]), "\n")
  fx <- warns(fixef(fit))
  print(fx$value)
}
